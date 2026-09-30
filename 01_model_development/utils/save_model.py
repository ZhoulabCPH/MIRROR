"""Persist model and optimizer states without serializing local storage paths."""

from pathlib import Path


def save_model(args, model, optimizer, current_epoch):
    """Save a checkpoint under ``out_dir/checkpoints`` and return its path."""
    import torch

    checkpoint_dir = Path(args.out_dir) / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = checkpoint_dir / f"checkpoint_{int(current_epoch):04d}.pt"
    state = {
        "net": model.state_dict(),
        "optimizer": optimizer.state_dict(),
        "epoch": int(current_epoch),
    }
    torch.save(state, checkpoint_path)
    return checkpoint_path
