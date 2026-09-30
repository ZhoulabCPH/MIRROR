"""Compute MIRROR prediction scores, numerical metrics, and export tables."""

import warnings

import numpy as np
import pandas as pd
import torch
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve


MODALITIES = ("Multimodal", "WSI", "RNASeq")
DATASET_ALIASES = (
    ("set_01", "partition_01"),
    ("set_02", "partition_02"),
    ("set_03", "partition_03"),
    ("set_04", "partition_04"),
    ("set_05", "partition_05"),
    ("set_06", "partition_06"),
)
EMPTY_RESULT = {
    "n": 0, "n_lab": 0, "n_event": 0, "n_pre0": 0, "n_pre1": 0,
    "acc": -1.0, "auc": -1.0, "sen": -1.0, "spe": -1.0,
    "c_index_hrd": -1.0, "p_logrank": -1.0, "hr": -1.0,
}


def youden_threshold(probabilities, labels):
    """Estimate a probability threshold from finite labeled observations."""
    probabilities, labels = np.asarray(probabilities), np.asarray(labels)
    valid = np.isfinite(probabilities) & np.isfinite(labels)
    if valid.sum() < 5 or len(np.unique(labels[valid])) < 2:
        return 0.5
    fpr, tpr, thresholds = roc_curve(labels[valid].astype(int), probabilities[valid])
    usable = np.isfinite(thresholds) & (thresholds >= 0) & (thresholds <= 1)
    if not usable.any():
        return 0.5
    indices = np.flatnonzero(usable)
    best_index = indices[int(np.argmax((tpr - fpr)[usable]))]
    return float(thresholds[best_index])


def compute_metrics(scores, probabilities, labels, times, events, threshold):
    """Compute classification and descriptive survival metrics at a fixed threshold.

    The concordance index treats higher HRD scores as longer predicted survival.
    A value of -1 denotes an unavailable metric.
    """
    scores, probabilities, labels, times, events = [
        np.asarray(values, dtype=float)
        for values in (scores, probabilities, labels, times, events)
    ]
    result = dict(EMPTY_RESULT)
    result["n"] = len(scores)
    finite_probability = np.isfinite(probabilities)
    prediction = probabilities[finite_probability] >= threshold
    result["n_pre0"] = int((~prediction).sum())
    result["n_pre1"] = int(prediction.sum())
    labeled = finite_probability & np.isfinite(labels)
    result["n_lab"] = int(labeled.sum())
    if labeled.any():
        truth = labels[labeled].astype(int)
        predicted = (probabilities[labeled] >= threshold).astype(int)
        result["acc"] = float(accuracy_score(truth, predicted))
        tn, fp, fn, tp = confusion_matrix(truth, predicted, labels=[0, 1]).ravel()
        result["sen"] = float(tp / (tp + fn)) if tp + fn else -1.0
        result["spe"] = float(tn / (tn + fp)) if tn + fp else -1.0
        if len(np.unique(truth)) == 2:
            result["auc"] = float(roc_auc_score(truth, probabilities[labeled]))

    survival = (
        np.isfinite(times) & (times > 0) & np.isin(events, [0, 1])
        & np.isfinite(scores) & finite_probability
    )
    result["n_event"] = int((events[survival] == 1).sum())
    if survival.sum() < 10 or result["n_event"] < 2:
        return result
    duration, observed = times[survival], events[survival]
    try:
        result["c_index_hrd"] = float(concordance_index(duration, scores[survival], observed))
    except (ValueError, ZeroDivisionError) as error:
        warnings.warn(f"Concordance is unavailable: {error}", RuntimeWarning)
    group = probabilities[survival] >= threshold
    if len(np.unique(group)) == 2:
        comparison = logrank_test(
            duration[group], duration[~group],
            event_observed_A=observed[group], event_observed_B=observed[~group],
        )
        result["p_logrank"] = float(comparison.p_value)
        try:
            frame = pd.DataFrame({"duration": duration, "event": observed, "group": group.astype(int)})
            fit = CoxPHFitter().fit(frame, duration_col="duration", event_col="event")
            result["hr"] = float(fit.hazard_ratios_["group"])
        except (ConvergenceError, ValueError, np.linalg.LinAlgError) as error:
            warnings.warn(f"Cox hazard ratio is unavailable: {error}", RuntimeWarning)
    return result


@torch.no_grad()
def collect_hrd_scores(model, dataset, device, modality="Multimodal"):
    """Collect per-sample scores with dropout disabled and gradients disabled."""
    if modality not in MODALITIES:
        raise ValueError(f"Unsupported modality: {modality}")
    model.eval()
    scores, probabilities, labels, times, events, identifiers = [], [], [], [], [], []
    for index in range(len(dataset)):
        sample = dataset[index]
        has_wsi = int(sample["mask_wsi"]) if modality != "RNASeq" else 0
        has_rna = int(sample["mask_SeqRNA"]) if modality != "WSI" else 0
        if not (has_wsi or has_rna):
            continue
        padding = sample.get("pad_mask")
        output = model(
            x_wsi=sample["WSI_feature"].unsqueeze(0).to(device),
            x_SeqRNA=sample["RNASeq_feature"].unsqueeze(0).to(device),
            mask_wsi=torch.tensor([has_wsi], device=device),
            mask_SeqRNA=torch.tensor([has_rna], device=device),
            pad_mask=None if padding is None else padding.unsqueeze(0).to(device),
            batch_size=1,
        )
        logit = output["HRD_score"][0]
        if not torch.isfinite(logit):
            raise FloatingPointError("Evaluation produces a non-finite HRD score.")
        scores.append(float(logit))
        probabilities.append(float(torch.sigmoid(logit)))
        labels.append(float(sample["Label"]))
        times.append(float(sample["OS"]))
        events.append(float(sample["OSState"]))
        identifiers.append(sample["patient_name"])
    return (
        np.asarray(scores), np.asarray(probabilities), np.asarray(labels),
        np.asarray(times), np.asarray(events), identifiers,
    )


def evaluate_checkpoint(model, datasets, device):
    """Compute per-modality scores and associated probability thresholds."""
    results, thresholds = {}, {}
    for modality in MODALITIES:
        scores = collect_hrd_scores(model, datasets["partition_01"], device, modality)
        threshold = youden_threshold(scores[1], scores[2])
        thresholds[modality] = threshold
        results[modality] = {
            "set_01": compute_metrics(*scores[:5], threshold),
        }
    return results, thresholds


def checkpoint_score(results):
    """Return the finite multimodal score used by the checkpoint comparator."""
    value = float(results["Multimodal"]["set_01"]["auc"])
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("Checkpoint scoring requires a finite AUC and both HRD classes.")
    return value


def export_final_evaluation(model, datasets, device, thresholds, output_dir):
    """Export predictions and metrics with generic dataset aliases."""
    predictions, metric_rows = [], []
    for modality in MODALITIES:
        threshold = float(thresholds[modality])
        for name, key in DATASET_ALIASES:
            collected = collect_hrd_scores(model, datasets[key], device, modality)
            scores, probabilities, labels, times, events, identifiers = collected
            metrics = compute_metrics(*collected[:5], threshold)
            metric_rows.append({"dataset": name, "modality": modality, "threshold": threshold, **metrics})
            predictions.append(pd.DataFrame({
                "sample_id": identifiers, "HRD_score": scores, "HRD_prob": probabilities,
                "PRE_HRD": (probabilities >= threshold).astype(int), "Label": labels,
                "OStime": times, "OS": events, "dataset": name, "modality": modality,
            }))
    pd.concat(predictions, ignore_index=True).to_csv(output_dir / "predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(output_dir / "evaluation_metrics.csv", index=False)
