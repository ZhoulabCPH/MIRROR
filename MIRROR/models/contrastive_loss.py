"""Instance, cluster, and overlap objectives for optional representation learning."""

import math

import torch
import torch.nn as nn


def _negative_pair_mask(size, device=None):
    """Exclude self-pairs and corresponding cross-view positive pairs."""
    mask = torch.ones((2 * size, 2 * size), dtype=torch.bool, device=device)
    mask.fill_diagonal_(False)
    indices = torch.arange(size, device=device)
    mask[indices, size + indices] = False
    mask[size + indices, indices] = False
    return mask


class InstanceLoss(nn.Module):
    """Contrast two views using corresponding rows as positive pairs.

    Embeddings are supplied in matching patient order. Their normalization is
    the caller's responsibility. The current input size supports partial batches.
    """

    def __init__(self, batch_size, temperature, device=None):
        super().__init__()
        if batch_size < 1 or temperature <= 0:
            raise ValueError("batch_size and temperature must be positive.")
        self.batch_size = batch_size
        self.temperature = temperature
        self.device = device
        self.register_buffer("mask", self.mask_correlated_samples(batch_size),
                             persistent=False)
        self.criterion = nn.CrossEntropyLoss(reduction="sum")

    def mask_correlated_samples(self, batch_size):
        """Return the negative-pair mask for a two-view batch."""
        return _negative_pair_mask(batch_size)

    def forward(self, z_i, z_j):
        """Return the mean instance contrastive loss for matched view embeddings."""
        if z_i.ndim != 2 or z_i.shape != z_j.shape or z_i.shape[0] == 0:
            raise ValueError("Two nonempty view matrices with identical shapes are required.")
        size = z_i.shape[0]
        count = 2 * size
        z = torch.cat((z_i, z_j), dim=0)
        similarity = z @ z.T / self.temperature
        positives = torch.cat((torch.diag(similarity, size),
                               torch.diag(similarity, -size))).reshape(count, 1)
        mask = (self.mask.to(z.device) if size == self.batch_size
                else _negative_pair_mask(size, z.device))
        negatives = similarity[mask].reshape(count, count - 2)
        logits = torch.cat((positives, negatives), dim=1)
        labels = torch.zeros(count, device=z.device, dtype=torch.long)
        return self.criterion(logits, labels) / count


class DiceLoss(nn.Module):
    """Measure one minus the soft Dice overlap of two equally shaped tensors."""

    def __init__(self, smooth=1e-6):
        super().__init__()
        if smooth <= 0:
            raise ValueError("smooth must be positive.")
        self.smooth = smooth

    def forward(self, y_pred, y_true):
        """Compute soft overlap across all entries, including noncontiguous tensors."""
        if y_pred.shape != y_true.shape:
            raise ValueError("Prediction and target shapes must match.")
        prediction, target = y_pred.reshape(-1), y_true.reshape(-1)
        intersection = torch.sum(prediction * target)
        dice = ((2.0 * intersection + self.smooth)
                / (prediction.sum() + target.sum() + self.smooth))
        return 1 - dice


class ClusterLoss(nn.Module):
    """Align cluster assignments across views while discouraging collapsed marginals."""

    def __init__(self, class_num, temperature, device=None):
        super().__init__()
        if class_num < 1 or temperature <= 0:
            raise ValueError("class_num and temperature must be positive.")
        self.class_num = class_num
        self.temperature = temperature
        self.device = device
        self.register_buffer("mask", self.mask_correlated_clusters(class_num),
                             persistent=False)
        self.criterion = nn.CrossEntropyLoss(reduction="sum")
        self.similarity_f = nn.CosineSimilarity(dim=2)

    def mask_correlated_clusters(self, class_num):
        """Exclude each cluster's self-pair and corresponding cluster in the other view."""
        return _negative_pair_mask(class_num)

    @staticmethod
    def _entropy_penalty(assignments):
        """Compute divergence from a uniform cluster marginal with finite zero terms."""
        totals = assignments.sum(dim=0)
        mass = totals.sum()
        if torch.any(assignments < 0) or not torch.isfinite(assignments).all() or mass <= 0:
            raise ValueError("Cluster assignments must be finite, nonnegative, and nonempty.")
        probabilities = totals / mass
        log_probabilities = probabilities.clamp_min(torch.finfo(probabilities.dtype).tiny).log()
        return math.log(probabilities.numel()) + (probabilities * log_probabilities).sum()

    def forward(self, c_i, c_j):
        """Return cluster contrastive loss plus marginal entropy penalties."""
        if (c_i.ndim != 2 or c_i.shape != c_j.shape
                or c_i.shape[1] != self.class_num or c_i.shape[0] == 0):
            raise ValueError("Matched assignment matrices must have class_num columns.")
        entropy = self._entropy_penalty(c_i) + self._entropy_penalty(c_j)
        count = 2 * self.class_num
        clusters = torch.cat((c_i.t(), c_j.t()), dim=0)
        similarity = self.similarity_f(clusters.unsqueeze(1), clusters.unsqueeze(0))
        similarity = similarity / self.temperature
        positives = torch.cat((torch.diag(similarity, self.class_num),
                               torch.diag(similarity, -self.class_num))).reshape(count, 1)
        negatives = similarity[self.mask.to(clusters.device)].reshape(count, count - 2)
        logits = torch.cat((positives, negatives), dim=1)
        labels = torch.zeros(count, device=clusters.device, dtype=torch.long)
        return self.criterion(logits, labels) / count + entropy
