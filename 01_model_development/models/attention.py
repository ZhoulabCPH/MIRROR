"""Reusable attention pooling and classifiers with optional confounder features.

Confounder feature files are caller-supplied NumPy arrays. The module loads
these embeddings and combines them with instance representations through
scaled attention.
"""

import math
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def _load_confounders(paths, feature_dim):
    """Load caller-selected feature files into one float tensor."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    features = [torch.from_numpy(np.load(path, allow_pickle=False))
                .reshape(-1, feature_dim).float() for path in paths]
    if not features:
        raise ValueError("At least one confounder feature file is required.")
    return torch.cat(features, dim=0)


class Attention2(nn.Module):
    """Generate independent sigmoid attention scores for each instance."""

    def __init__(self, L=512, D=128, K=1):
        super().__init__()
        self.L, self.D, self.K = L, D, K
        self.attention = nn.Sequential(
            nn.Linear(L, D), nn.Tanh(), nn.Linear(D, K)
        )

    def forward(self, x, isNorm=True):
        """Return ``(K, instances)`` scores, optionally bounded by a sigmoid."""
        scores = self.attention(x).transpose(1, 0)
        return scores.sigmoid() if isNorm else scores


class Attention_Gated(nn.Module):
    """Generate gated attention with parameter sharing across pooling heads.

    All heads share the same feature projections and score layer. Consequently,
    each head produces the same scores for a given input. This is a shared-head
    pooling module, not independent multihead self-attention.
    """

    def __init__(self, L=512, D=128, K=1, num_heads=8):
        super().__init__()
        if num_heads < 1:
            raise ValueError("num_heads must be positive.")
        self.L, self.D, self.K = L, D, K
        self.num_heads = num_heads
        self.attention_V = nn.Sequential(nn.Linear(L, D), nn.Tanh())
        self.attention_U = nn.Sequential(nn.Linear(L, D), nn.Sigmoid())
        self.attention_weights = nn.Linear(D, K)
        # Shared references retain the intended parameter-sharing structure.
        self.attention_V_heads = nn.ModuleList([self.attention_V for _ in range(num_heads)])
        self.attention_U_heads = nn.ModuleList([self.attention_U for _ in range(num_heads)])
        self.attention_weights_heads = nn.ModuleList(
            [self.attention_weights for _ in range(num_heads)]
        )

    def forward(self, x, isNorm=True):
        """Return per-head attention and the mean weight for each patch."""
        scores = []
        for head in range(self.num_heads):
            gated = self.attention_V_heads[head](x) * self.attention_U_heads[head](x)
            scores.append(self.attention_weights_heads[head](gated).transpose(1, 0))
        attention = torch.cat(scores, dim=0)
        if isNorm:
            attention = F.softmax(attention, dim=1)
        return attention, attention.mean(dim=0)


class Classifier_1fc(nn.Module):
    """Classify embeddings with optional attention over confounder features."""

    def __init__(self, n_channels, n_classes, droprate=0.0, confounder_path=False):
        super().__init__()
        self.confounder_path = confounder_path
        self.droprate = droprate
        if droprate != 0.0:
            self.dropout = nn.Dropout(p=droprate)
        if confounder_path:
            self.register_buffer("confounder_feat", _load_confounders(confounder_path, n_channels))
            self.W_q = nn.Linear(n_channels, 128)
            self.W_k = nn.Linear(n_channels, 128)
            self.fc = nn.Linear(n_channels * 2, n_classes)
        else:
            self.fc = nn.Linear(n_channels, n_classes)

    def forward(self, x):
        """Return logits, classifier embeddings, and optional confounder attention."""
        if self.droprate != 0.0:
            x = self.dropout(x)
        if self.confounder_path:
            query = self.W_q(x)
            keys = self.W_k(self.confounder_feat)
            scores = keys @ query.transpose(0, 1) / math.sqrt(keys.shape[1])
            attention = F.softmax(scores, dim=0)
            confounders = attention.transpose(0, 1) @ self.confounder_feat
            embedding = torch.cat((x, confounders), dim=1)
            return self.fc(embedding), embedding, attention
        return self.fc(x), x, None


class Attention_with_Classifier(nn.Module):
    """Pool an instance bag and classify each pooled representation.

    ``args`` is accepted for caller compatibility; no configuration fields are
    read from it. Optional confounder arrays contain features with dimension
    ``L``.
    """

    def __init__(self, args, L=512, D=128, K=1, num_cls=2,
                 droprate=0, confounder_path=False):
        super().__init__()
        self.attention = Attention_Gated(L, D, K)
        self.confounder_path = confounder_path
        if confounder_path:
            self.register_buffer("confounder_feat", _load_confounders(confounder_path, L))
            self.W_q = nn.Linear(L, 128)
            self.W_k = nn.Linear(L, 128)
            self.classifier = nn.Linear(L * 2, num_cls)
            self.dropout = nn.Dropout(0.5)
        else:
            self.classifier = Classifier_1fc(L, num_cls, droprate)

    def forward(self, x):
        """Return logits, pooled embeddings, and the attention used for classification."""
        attention, _ = self.attention(x)
        embedding = attention @ x
        if self.confounder_path:
            query = self.W_q(embedding)
            keys = self.W_k(self.confounder_feat)
            scores = keys @ query.transpose(0, 1) / math.sqrt(keys.shape[1])
            confounder_attention = F.softmax(scores, dim=0)
            confounders = confounder_attention.transpose(0, 1) @ self.confounder_feat
            embedding = torch.cat((embedding, confounders), dim=1)
            return self.classifier(embedding), embedding, confounder_attention
        logits, _, _ = self.classifier(embedding)
        return logits, embedding, attention
