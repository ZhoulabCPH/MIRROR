"""MIRROR encoders, missing-modality fusion, and objective components."""

from typing import Dict

import torch
import torch.nn as nn
import torch.nn.functional as F


class GatedAttentionMIL(nn.Module):
    """Pool patch features into slide embeddings with gated attention.

    Inputs have shape ``(batch, patches, in_dim)``. The optional patch mask
    marks observed patches with one and padding with zero. Padded patches
    receive zero attention, including bags that contain only padding.
    """

    def __init__(self, in_dim=1024, d=256, att_dim=128, dropout=0.25):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(in_dim, d), nn.LayerNorm(d), nn.GELU(), nn.Dropout(dropout),
        )
        self.att_V = nn.Sequential(nn.Linear(d, att_dim), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(d, att_dim), nn.Sigmoid())
        self.att_w = nn.Linear(att_dim, 1)
        self.norm = nn.LayerNorm(d)

    def forward(self, x, pad_mask=None):
        """Return slide embeddings and one attention weight per patch."""
        if x.ndim != 3 or x.shape[1] == 0:
            raise ValueError("WSI features require a nonempty (batch, patches, features) tensor.")
        h = self.proj(x)
        scores = self.att_w(self.att_V(h) * self.att_U(h))
        valid = None
        if pad_mask is not None:
            valid = torch.as_tensor(pad_mask, device=x.device)
            if valid.shape != x.shape[:2]:
                raise ValueError("pad_mask must match the WSI batch and patch dimensions.")
            if not torch.all((valid == 0) | (valid == 1)):
                raise ValueError("pad_mask contains values other than zero or one.")
            valid = valid.bool().unsqueeze(-1)
            scores = scores.masked_fill(~valid, torch.finfo(scores.dtype).min)
        weights = torch.softmax(scores, dim=1)
        if valid is not None:
            weights = weights * valid
            weights = weights / weights.sum(dim=1, keepdim=True).clamp_min(
                torch.finfo(weights.dtype).tiny
            )
        pooled = torch.sum(weights * h, dim=1)
        return self.norm(pooled), weights.squeeze(-1)


class RNAEncoder(nn.Module):
    """Map a fixed gene-expression vector to a normalized RNA embedding."""

    def __init__(self, in_dim=12042, d=256, hidden=512, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(hidden, d), nn.LayerNorm(d), nn.GELU(),
        )

    def forward(self, x):
        """Encode one gene-expression vector per patient."""
        return self.net(x)


class MIRROR(nn.Module):
    """Predict HRD from whole-slide and RNA features with missing-modality tokens.

    ``HRD_score`` is the prediction logit; larger values indicate greater HRD
    likelihood. Auxiliary logits and normalized embeddings support supervised
    unimodal objectives and cross-modal alignment. Attention
    weights support patch-level interpretation when slide features are present.
    """

    def __init__(self, d_wsi: int = 1024, d_rna: int = 12042,
                 d_model: int = 256, d_fusion: int = 128,
                 proj_dim: int = 128, dropout: float = 0.25):
        super().__init__()
        self.d_model = d_model
        self.wsi_mil = GatedAttentionMIL(d_wsi, d_model, att_dim=128, dropout=dropout)
        self.rna_enc = RNAEncoder(d_rna, d_model, hidden=512, dropout=0.3)

        # Learned tokens represent absent modalities before fusion.
        self.token_wsi = nn.Parameter(torch.zeros(1, d_model))
        self.token_rna = nn.Parameter(torch.zeros(1, d_model))
        nn.init.normal_(self.token_wsi, std=0.02)
        nn.init.normal_(self.token_rna, std=0.02)

        self.cross_w2r = nn.MultiheadAttention(
            d_model, num_heads=4, dropout=dropout, batch_first=True
        )
        self.cross_r2w = nn.MultiheadAttention(
            d_model, num_heads=4, dropout=dropout, batch_first=True
        )
        self.ln_w = nn.LayerNorm(d_model)
        self.ln_r = nn.LayerNorm(d_model)
        self.fusion = nn.Sequential(
            nn.LayerNorm(2 * d_model),
            nn.Linear(2 * d_model, d_fusion), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(d_fusion, d_fusion), nn.GELU(),
        )
        self.hrd_head = nn.Sequential(nn.Dropout(dropout), nn.Linear(d_fusion, 1))

        # Unimodal heads and projection layers supply auxiliary outputs.
        self.hrd_head_wsi = nn.Linear(d_model, 1)
        self.hrd_head_rna = nn.Linear(d_model, 1)
        self.proj_wsi = nn.Sequential(nn.Linear(d_model, proj_dim), nn.GELU(),
                                      nn.Linear(proj_dim, proj_dim))
        self.proj_rna = nn.Sequential(nn.Linear(d_model, proj_dim), nn.GELU(),
                                      nn.Linear(proj_dim, proj_dim))

    @staticmethod
    def _as_batch_tensor(x, device):
        """Convert stacked tensors or equal-sized sample lists to model inputs."""
        if x is None:
            return None
        if isinstance(x, (list, tuple)):
            x = torch.stack([torch.as_tensor(sample) for sample in x], dim=0)
        return torch.as_tensor(x, device=device, dtype=torch.float32)

    @staticmethod
    def _modality_mask(mask, features, size, device, name):
        """Validate binary availability and infer absence from a missing input."""
        if mask is None:
            return torch.full((size,), float(features is not None), device=device)
        mask = torch.as_tensor(mask, device=device, dtype=torch.float32).reshape(-1)
        if mask.numel() != size or not torch.all((mask == 0) | (mask == 1)):
            raise ValueError(f"{name} must contain one binary value per patient.")
        if features is None and mask.any():
            raise ValueError(f"{name} marks an absent input as available.")
        return mask

    def forward(self, x_wsi=None, x_SeqRNA=None, mask_wsi=None,
                mask_SeqRNA=None, pad_mask=None,
                batch_size=None) -> Dict[str, torch.Tensor]:
        """Predict a batch using available modalities and optional patch masks.

        WSI features have shape ``(B, N, d_wsi)`` and RNA features have shape
        ``(B, d_rna)``. Modality masks contain ``B`` binary availability values.
        ``pad_mask`` has shape ``(B, N)`` and identifies real patches. A single
        unbatched sample is accepted for either modality. ``batch_size``, when
        supplied, validates the inferred number of patients.
        """
        device = next(self.parameters()).device
        x_wsi = self._as_batch_tensor(x_wsi, device)
        x_SeqRNA = self._as_batch_tensor(x_SeqRNA, device)
        if x_wsi is not None and x_wsi.ndim == 2:
            x_wsi = x_wsi.unsqueeze(0)
        if x_SeqRNA is not None and x_SeqRNA.ndim == 1:
            x_SeqRNA = x_SeqRNA.unsqueeze(0)
        if x_wsi is None and x_SeqRNA is None:
            raise ValueError("At least one of x_wsi and x_SeqRNA must be provided.")
        if x_wsi is not None and x_wsi.ndim != 3:
            raise ValueError("WSI features must have shape (batch, patches, features).")
        if x_SeqRNA is not None and x_SeqRNA.ndim != 2:
            raise ValueError("RNA features must have shape (batch, genes).")
        size = x_wsi.shape[0] if x_wsi is not None else x_SeqRNA.shape[0]
        if x_SeqRNA is not None and x_SeqRNA.shape[0] != size:
            raise ValueError("WSI and RNA features must have matching batch sizes.")
        if batch_size is not None and int(batch_size) != size:
            raise ValueError("batch_size does not match the supplied feature tensors.")
        m_w = self._modality_mask(mask_wsi, x_wsi, size, device, "mask_wsi")
        m_r = self._modality_mask(mask_SeqRNA, x_SeqRNA, size, device, "mask_SeqRNA")

        if x_wsi is not None:
            h_w, attention = self.wsi_mil(x_wsi, pad_mask)
            if pad_mask is not None:
                # Bags containing only padding carry no observed slide modality.
                m_w = m_w * torch.as_tensor(pad_mask, device=device).bool().any(dim=1)
            attention = attention * m_w.unsqueeze(1)
        else:
            h_w = torch.zeros(size, self.d_model, device=device)
            attention = None
        h_r = (self.rna_enc(x_SeqRNA) if x_SeqRNA is not None
               else torch.zeros(size, self.d_model, device=device))

        h_w = torch.where(m_w.bool().unsqueeze(1), h_w, self.token_wsi.expand(size, -1))
        h_r = torch.where(m_r.bool().unsqueeze(1), h_r, self.token_rna.expand(size, -1))

        # Each direction attends to the single embedding of the other modality.
        qw, qr = h_w.unsqueeze(1), h_r.unsqueeze(1)
        aw, _ = self.cross_w2r(qw, qr, qr)
        ar, _ = self.cross_r2w(qr, qw, qw)
        f_w = self.ln_w(h_w + aw.squeeze(1))
        f_r = self.ln_r(h_r + ar.squeeze(1))
        fused = self.fusion(torch.cat([f_w, f_r], dim=1))
        outputs = {
            "HRD_score": self.hrd_head(fused).reshape(-1),
            "hrd_wsi": self.hrd_head_wsi(h_w).reshape(-1),
            "hrd_rna": self.hrd_head_rna(h_r).reshape(-1),
            "z_wsi": F.normalize(self.proj_wsi(h_w), dim=1),
            "z_rna": F.normalize(self.proj_rna(h_r), dim=1),
            "mask_wsi": m_w,
            "mask_SeqRNA": m_r,
        }
        if attention is not None:
            outputs["attention"] = attention
        return outputs


def masked_bce(logit, label, valid, pos_weight=None):
    """Compute binary cross-entropy for valid labels and preserve zero gradients."""
    valid = torch.as_tensor(valid, device=logit.device).bool()
    if not valid.any():
        return logit.sum() * 0.0
    return F.binary_cross_entropy_with_logits(
        logit[valid], label[valid].float(), pos_weight=pos_weight
    )


def info_nce(z_a, z_b, valid, temperature=0.1):
    """Align paired observed modalities with symmetric instance discrimination."""
    if temperature <= 0:
        raise ValueError("temperature must be positive.")
    valid = torch.as_tensor(valid, device=z_a.device).bool()
    count = int(valid.sum().item())
    if count < 2:
        return (z_a.sum() + z_b.sum()) * 0.0
    logits = z_a[valid] @ z_b[valid].t() / temperature
    target = torch.arange(count, device=z_a.device)
    return 0.5 * (F.cross_entropy(logits, target)
                  + F.cross_entropy(logits.t(), target))
