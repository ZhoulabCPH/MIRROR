"""Optimize MIRROR, manage checkpoints, and export model outputs."""

import copy
import hashlib
import json
import random
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import torch
from torch.optim.lr_scheduler import CosineAnnealingLR

from datasets.data_loaders import build_data_loaders
from evaluation import evaluate_checkpoint, export_final_evaluation, checkpoint_score
from models.mirror import MIRROR, info_nce, masked_bce
from utils.config import parse_config
from utils.metric_log import MetricLogger


PIPELINE_VERSION = "mirror_pipeline_v2"
SCORING_VERSION = "mirror_score_v2"


def set_seed(seed):
    """Seed the Python, NumPy, and PyTorch generators for repeatable sampling."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def move_batch_to_device(batch, device):
    """Move tensor fields to the selected device and retain sample metadata."""
    return {key: value.to(device) if torch.is_tensor(value) else value for key, value in batch.items()}


def compute_loss(outputs, batch, args):
    """Combine supervised HRD loss, unimodal supervision, and paired alignment."""
    label = batch["Label"].reshape(-1).float()
    if not torch.isfinite(label).all() or not ((label == 0) | (label == 1)).all():
        raise ValueError("Every observation requires a binary HRD label.")
    has_wsi, has_rna = outputs["mask_wsi"] > 0.5, outputs["mask_SeqRNA"] > 0.5
    if not (has_wsi | has_rna).all():
        raise ValueError("Every observation requires at least one available modality.")
    classification = masked_bce(outputs["HRD_score"], label, torch.ones_like(label))
    unimodal = (
        masked_bce(outputs["hrd_wsi"], label, has_wsi)
        + masked_bce(outputs["hrd_rna"], label, has_rna)
    )
    alignment = info_nce(outputs["z_wsi"], outputs["z_rna"], has_wsi & has_rna)
    total = args.w_cls * classification + args.w_uni * unimodal + args.w_align * alignment
    return total, {
        "total": float(total.detach()), "cls": float(classification.detach()),
        "uni": float(unimodal.detach()), "align": float(alignment.detach()),
    }


def assert_input_partition(dataset):
    """Check the required metadata conditions before processing batches."""
    frame = dataset.df
    if frame.empty or not ((frame["Strategy"] == "PARTITION_01") & (frame["cohort"] == "SOURCE_A")).all():
        raise ValueError("The input partition is empty or contains incompatible metadata.")


def partition_fingerprint(datasets):
    """Hash sample membership and labels without exporting identifiers."""
    rows = []
    for key in ("partition_01",):
        frame = datasets[key].df
        for _, row in frame.iterrows():
            rows.append((key, str(row["sample_id"]), str(row["HRD_status"])))
    return hashlib.sha256(json.dumps(sorted(rows), separators=(",", ":")).encode("utf-8")).hexdigest()


def load_checkpoint(path, device, fingerprint):
    """Load a checkpoint with compatible metadata and sample fingerprints."""
    checkpoint = torch.load(path, map_location=device, weights_only=True)
    required = {"state_dict", "optimizer", "scheduler", "best", "rng", "epoch"}
    if not isinstance(checkpoint, dict) or not required.issubset(checkpoint):
        raise ValueError("Continuation requires a best_model.pth or last_model.pth resumable checkpoint.")
    if checkpoint.get("pipeline_version") != PIPELINE_VERSION:
        raise ValueError("The checkpoint pipeline version is incompatible.")
    if checkpoint.get("scoring_version") != SCORING_VERSION:
        raise ValueError("The checkpoint scoring version is incompatible.")
    if checkpoint.get("partition_fingerprint") != fingerprint:
        raise ValueError("The checkpoint sample membership or labels do not match this run.")
    return checkpoint


def checkpoint_payload(model, optimizer, scheduler, loader, epoch, fingerprint, best):
    """Capture optimizer state and the selected model without local input paths."""
    numpy_rng = np.random.get_state()
    return {
        "model_name": "MIRROR", "epoch": epoch,
        "state_dict": model.state_dict(), "optimizer": optimizer.state_dict(),
        "scheduler": scheduler.state_dict(), "partition_fingerprint": fingerprint,
        "pipeline_version": PIPELINE_VERSION, "scoring_version": SCORING_VERSION,
        "best": best,
        "rng": {
            "python": random.getstate(),
            "numpy": (numpy_rng[0], numpy_rng[1].tolist(), int(numpy_rng[2]), int(numpy_rng[3]), float(numpy_rng[4])),
            "torch": torch.get_rng_state(),
            "cuda": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else [],
            "loader": loader.generator.get_state(),
        },
    }


def restore_rng(checkpoint, loader):
    """Restore sampling generators from the checkpoint at an epoch boundary."""
    state = checkpoint["rng"]
    random.setstate(state["python"])
    numpy_rng = state["numpy"]
    np.random.set_state((numpy_rng[0], np.asarray(numpy_rng[1], dtype=np.uint32), *numpy_rng[2:]))
    torch.set_rng_state(state["torch"].cpu())
    if torch.cuda.is_available() and state["cuda"]:
        if len(state["cuda"]) != torch.cuda.device_count():
            raise ValueError("Checkpoint continuation requires the same CUDA device count.")
        torch.cuda.set_rng_state_all([value.cpu() for value in state["cuda"]])
    loader.generator.set_state(state["loader"].cpu())


def optimize_epoch(model, loader, optimizer, device, args):
    """Update model weights and accumulate sample-weighted loss components."""
    assert_input_partition(loader.dataset)
    model.train()
    totals = {"total": 0.0, "cls": 0.0, "uni": 0.0, "align": 0.0}
    count = 0
    for batch in loader:
        batch = move_batch_to_device(batch, device)
        optimizer.zero_grad(set_to_none=True)
        outputs = model(
            x_wsi=batch["WSI_feature"], x_SeqRNA=batch["RNASeq_feature"],
            mask_wsi=batch["mask_wsi"], mask_SeqRNA=batch["mask_SeqRNA"],
            pad_mask=batch.get("pad_mask"), batch_size=len(batch["Label"]),
        )
        loss, components = compute_loss(outputs, batch, args)
        if not torch.isfinite(loss):
            raise FloatingPointError("The objective is non-finite.")
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0, error_if_nonfinite=True)
        optimizer.step()
        for name, value in components.items():
            totals[name] += value * len(batch["Label"])
        count += len(batch["Label"])
    if count == 0:
        raise ValueError("The data loader does not yield any observations.")
    return {name: value / count for name, value in totals.items()}


def main(argv=None):
    """Create a MIRROR run and persist model state and prediction outputs."""
    args = parse_config(argv, description="Optimize MIRROR and export model outputs.")
    set_seed(int(args.seed))
    epochs, evaluation_interval = int(args.Epoch), int(args.eval_every)
    if epochs < 1 or evaluation_interval < 1:
        raise ValueError("Epoch and eval_every must be positive integers.")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    batch_loader, datasets = build_data_loaders(args)
    assert_input_partition(batch_loader.dataset)
    score_labels = datasets["partition_01"].df["HRD_status"].dropna().unique()
    if set(score_labels) != {0, 1}:
        raise ValueError("The scoring input requires both HRD classes.")
    fingerprint = partition_fingerprint(datasets)
    model = MIRROR(d_wsi=int(args.D_WSI), d_rna=int(args.RNAseq)).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(args.start_lr), weight_decay=float(args.weight_decay))
    scheduler = CosineAnnealingLR(optimizer, T_max=epochs, eta_min=1e-6)
    best = {"auc": -1.0, "epoch": 0, "thresholds": {}, "state_dict": None}
    start_epoch = 0
    if args.initial_checkpoint is not None:
        checkpoint = load_checkpoint(args.initial_checkpoint, device, fingerprint)
        if int(checkpoint["scheduler"]["T_max"]) != epochs:
            raise ValueError("Checkpoint continuation requires the same configured Epoch horizon.")
        model.load_state_dict(checkpoint["state_dict"], strict=True)
        optimizer.load_state_dict(checkpoint["optimizer"])
        scheduler.load_state_dict(checkpoint["scheduler"])
        best = checkpoint["best"]
        start_epoch = int(checkpoint["epoch"])
        if start_epoch >= epochs:
            raise ValueError("Epoch must exceed the completed checkpoint epoch.")
        restore_rng(checkpoint, batch_loader)

    output_dir = Path(args.out_dir) / "MIRROR" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output_dir.mkdir(parents=True, exist_ok=False)
    logger = MetricLogger(output_dir / "run_metrics.docx")
    logger.log(f"MIRROR | device={device.type} | observations={len(batch_loader.dataset)}")

    for epoch in range(start_epoch, epochs):
        started = time.perf_counter()
        learning_rate = optimizer.param_groups[0]["lr"]
        losses = optimize_epoch(model, batch_loader, optimizer, device, args)
        scheduler.step()
        logger.log(
            f"Epoch {epoch + 1}/{epochs} | LR={learning_rate:.6g} | "
            f"loss={losses['total']:.4f}, cls={losses['cls']:.4f}, "
            f"uni={losses['uni']:.4f}, align={losses['align']:.4f} | "
            f"duration={time.perf_counter() - started:.1f}s"
        )
        if (epoch + 1) % evaluation_interval == 0 or epoch + 1 == epochs:
            results, thresholds = evaluate_checkpoint(model, datasets, device)
            candidate_auc = checkpoint_score(results)
            logger.log_results(results, thresholds)
            if candidate_auc > best["auc"]:
                best = {
                    "auc": candidate_auc, "epoch": epoch + 1, "thresholds": thresholds,
                    "state_dict": {name: value.detach().cpu().clone() for name, value in model.state_dict().items()},
                }
                torch.save(
                    checkpoint_payload(model, optimizer, scheduler, batch_loader, epoch + 1, fingerprint, copy.deepcopy(best)),
                    output_dir / "best_model.pth",
                )
                logger.log(f"Selected checkpoint: epoch={best['epoch']}, AUC={best['auc']:.4f}")
            torch.save(
                checkpoint_payload(model, optimizer, scheduler, batch_loader, epoch + 1, fingerprint, best),
                output_dir / "last_model.pth",
            )
            logger.save()

    if best["state_dict"] is None:
        raise RuntimeError("The run does not produce a valid selected checkpoint.")
    model.load_state_dict(best["state_dict"], strict=True)
    model.eval()
    model.requires_grad_(False)
    # Persist the selected parameters and associated numerical thresholds.
    final_checkpoint = {
        "model_name": "MIRROR", "epoch": best["epoch"], "state_dict": best["state_dict"],
        "pipeline_version": PIPELINE_VERSION, "scoring_version": SCORING_VERSION,
        "partition_fingerprint": fingerprint, "selection_auc": best["auc"],
        "hrd_threshold": best["thresholds"]["Multimodal"], "thresholds": best["thresholds"],
        "dimensions": {"d_wsi": int(args.D_WSI), "d_rna": int(args.RNAseq)},
    }
    torch.save(final_checkpoint, output_dir / "mirror_model.pth")
    with (output_dir / "thresholds.json").open("w", encoding="utf-8") as stream:
        json.dump(best["thresholds"], stream, indent=2)
    export_final_evaluation(model, datasets, device, best["thresholds"], output_dir)
    logger.log(f"Selected MIRROR epoch={best['epoch']}; final predictions and evaluation metrics are saved.")
    logger.save()


if __name__ == "__main__":
    main()
