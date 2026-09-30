"""Construct MIRROR data loaders and named dataset objects."""

import random

import numpy as np
import torch
from torch.utils.data import DataLoader

from .multimodal_dataset import (
    MirrorDataset,
    MirrorRecordDataset,
    collate_multimodal,
    load_cohorts,
)


def _seed_worker(worker_id):
    """Seed Python and NumPy from the worker seed assigned by PyTorch."""
    worker_seed = torch.initial_seed() % (2 ** 32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def build_data_loaders(args):
    """Return the configured loader and named dataset objects."""
    cohorts = load_cohorts(args)
    partitions = cohorts["partitions"]
    partition_01 = MirrorDataset(
        partitions["partition_01"], cohorts["rna"],
        cohorts["feature_directories"]["SOURCE_A"], args,
    )
    batch_size = int(args.batch_size)
    workers = int(args.workers)
    if batch_size <= 0 or workers < 0:
        raise ValueError("Batch size must be positive and worker count must be nonnegative.")
    generator = torch.Generator()
    generator.manual_seed(int(getattr(args, "seed", 42)))
    batch_loader = DataLoader(
        partition_01,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        num_workers=workers,
        pin_memory=False,
        collate_fn=collate_multimodal,
        generator=generator,
        worker_init_fn=_seed_worker,
    )
    dataset_bank = {
        name: MirrorRecordDataset(
            clinical, cohorts["rna"],
            cohorts["feature_directories"][cohorts["partition_cohorts"][name]], args,
        )
        for name, clinical in partitions.items()
    }
    args.dim = int(args.D_WSI)
    print(f"[Data] {len(partition_01)} samples, {len(batch_loader)} batches per epoch.")
    return batch_loader, dataset_bank
