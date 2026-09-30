"""Validate clinical metadata and load MIRROR multimodal features."""

from hashlib import sha256
import os
from pathlib import Path
import re
from tempfile import NamedTemporaryFile

import h5py
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset


SOURCE_PARTITIONS = {
    "SOURCE_A": frozenset({"PARTITION_01", "PARTITION_02", "PARTITION_03"}),
    "SOURCE_B": frozenset({"PARTITION_04"}),
    "SOURCE_C": frozenset({"PARTITION_05"}),
    "SOURCE_D": frozenset({"PARTITION_06"}),
}
PARTITION_KEYS = tuple(f"PARTITION_{index:02d}" for index in range(1, 7))
CLINICAL_COLUMNS = {
    "sample_id", "cohort", "IS_SeqRNA", "IS_FFPE", "HRD_status",
    "OS", "OStime", "DFS", "DFStime",
}
_WSI_FEATURE_CACHE = {}
_CACHE_VERSION = 1


def _identifier(value):
    """Normalize an identifier for case-insensitive overlap comparisons."""
    return str(value).strip().casefold()


def _participant_id(sample_id, participant_id_pattern=None):
    """Extract a participant identifier with a configured capture pattern."""
    if not participant_id_pattern:
        return None
    try:
        expression = re.compile(participant_id_pattern, re.I)
    except (re.error, TypeError):
        raise ValueError("The participant identifier pattern is invalid.") from None
    if expression.groups < 1:
        raise ValueError("The participant identifier pattern requires a capture group.")
    match = expression.match(str(sample_id))
    return _identifier(match.group(1)) if match and match.group(1) else None


def source_labels_from_args(args):
    """Read configurable source labels using generic public defaults."""
    return {
        source: getattr(args, f"{source.lower()}_label", source)
        for source in SOURCE_PARTITIONS
    }


def partition_labels_from_args(args):
    """Read configurable partition labels using generic public defaults."""
    return {
        partition: getattr(args, f"{partition.lower()}_label", partition)
        for partition in PARTITION_KEYS
    }


def _partition_aliases(partition_labels=None):
    """Resolve configured partition values to neutral internal keys."""
    labels = dict(partition_labels or {})
    if set(labels).difference(PARTITION_KEYS):
        raise ValueError("Partition label settings contain an unsupported key.")
    aliases = {_identifier(partition): partition for partition in PARTITION_KEYS}
    for partition in PARTITION_KEYS:
        label = labels.get(partition, partition)
        if label is None or not str(label).strip():
            raise ValueError("Partition labels must be nonempty.")
        alias = _identifier(label)
        if alias in aliases and aliases[alias] != partition:
            raise ValueError("Partition labels must identify distinct partitions.")
        aliases[alias] = partition
    return aliases


def _source_aliases(source_labels=None):
    """Build an unambiguous mapping from configured labels to source keys."""
    labels = dict(source_labels or {})
    if set(labels).difference(SOURCE_PARTITIONS):
        raise ValueError("Source label settings contain an unsupported source key.")
    aliases = {_identifier(source): source for source in SOURCE_PARTITIONS}
    for source in SOURCE_PARTITIONS:
        label = labels.get(source, source)
        if label is None or not str(label).strip():
            raise ValueError("Source labels must be nonempty.")
        alias = _identifier(label)
        if alias in aliases and aliases[alias] != source:
            raise ValueError("Source labels must identify distinct sources.")
        aliases[alias] = source
    return aliases


def _identity_columns(frame, patient_id_column=None):
    """Locate explicitly named patient identifiers in clinical metadata."""
    if patient_id_column:
        if patient_id_column not in frame.columns:
            raise ValueError("The configured patient identifier column is missing.")
        identifiers = frame[patient_id_column]
        if identifiers.isna().any() or identifiers.astype(str).str.strip().eq("").any():
            raise ValueError("The configured patient identifier column contains missing identifiers.")
        return [patient_id_column]
    return [name for name in ("patient_id", "case_id") if name in frame.columns]


def patient_group_keys(frame, patient_id_column=None, participant_id_pattern=None):
    """Resolve grouping keys from configured participant identifiers."""
    columns = _identity_columns(frame, patient_id_column)
    keys = []
    for _, row in frame.iterrows():
        key = _participant_id(row["sample_id"], participant_id_pattern)
        if key is None:
            key = next(
                (_identifier(row[column]) for column in columns
                 if pd.notna(row[column]) and str(row[column]).strip()),
                _identifier(row["sample_id"]),
            )
        keys.append(key)
    return pd.Series(keys, index=frame.index, dtype="string")


def validate_cohort_table(frame, patient_id_column=None, require_strategy=True,
                          source_labels=None, participant_id_pattern=None,
                          partition_labels=None):
    """Validate identifiers, modality metadata, and partition assignments."""
    required = CLINICAL_COLUMNS | ({"Strategy"} if require_strategy else set())
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Clinical metadata is missing required columns: {', '.join(missing)}.")
    frame = frame.copy().reset_index(drop=True)
    if frame.empty:
        raise ValueError("Clinical metadata contains no samples.")
    if frame["sample_id"].isna().any():
        raise ValueError("Clinical sample identifiers must not be missing.")
    frame["sample_id"] = frame["sample_id"].astype(str).str.strip()
    invalid_id = frame["sample_id"].map(
        lambda value: not value or value in {".", ".."}
        or any(character in value for character in '/\\:\x00')
    )
    if invalid_id.any():
        raise ValueError("Sample identifiers must be nonempty filename stems, without path separators.")
    if frame["sample_id"].map(_identifier).duplicated().any():
        raise ValueError("Duplicate sample identifiers occur in clinical metadata.")
    aliases = _source_aliases(source_labels)
    frame["cohort"] = frame["cohort"].map(
        lambda value: aliases.get(_identifier(value)) if pd.notna(value) else None
    )
    if (~frame["cohort"].isin(SOURCE_PARTITIONS)).any():
        raise ValueError("Clinical metadata contains a missing or unsupported cohort.")
    for column in ("IS_SeqRNA", "IS_FFPE"):
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if (~numeric.isin([0, 1])).any():
            raise ValueError(f"{column} contains values other than 0 or 1.")
        frame[column] = numeric.astype(np.int64)
    if ((frame["IS_SeqRNA"] == 0) & (frame["IS_FFPE"] == 0)).any():
        raise ValueError("Each clinical record must have at least one available modality.")
    for column in ("HRD_status", "OS", "DFS"):
        numeric = pd.to_numeric(frame[column], errors="coerce")
        invalid = frame[column].notna() & ~numeric.isin([0, 1])
        if invalid.any():
            raise ValueError(f"{column} must contain binary values or missing observations.")
        frame[column] = numeric.astype(float)
    for column in ("OStime", "DFStime"):
        numeric = pd.to_numeric(frame[column], errors="coerce")
        invalid = frame[column].notna() & (~np.isfinite(numeric) | (numeric < 0))
        if invalid.any():
            raise ValueError(f"{column} must contain finite nonnegative values or missing observations.")
        frame[column] = numeric.astype(float)
    if require_strategy:
        strategies = _partition_aliases(partition_labels)
        frame["Strategy"] = frame["Strategy"].map(
            lambda value: strategies.get(_identifier(value)) if pd.notna(value) else None
        )
        for cohort, allowed in SOURCE_PARTITIONS.items():
            if (~frame.loc[frame["cohort"] == cohort, "Strategy"].isin(allowed)).any():
                raise ValueError(f"Cohort/Strategy mismatch: {cohort} contains an invalid partition assignment.")
    # Shared identity tokens detect aliases across identifier columns.
    patient_columns = _identity_columns(frame, patient_id_column)
    assignments = {}
    for _, row in frame.iterrows():
        tokens = {_identifier(row["sample_id"])}
        participant = _participant_id(row["sample_id"], participant_id_pattern)
        if participant:
            tokens.add(participant)
        for column in patient_columns:
            if pd.notna(row[column]) and str(row[column]).strip():
                tokens.add(_identifier(row[column]))
        partition = (row["cohort"], row["Strategy"] if require_strategy else None)
        for token in tokens:
            previous = assignments.setdefault(token, partition)
            if previous != partition:
                raise ValueError("Patient overlap occurs across cohort or partition boundaries.")
    return frame


def _read_csv(path, description):
    """Read a table without including private paths in error messages."""
    if path is None or not str(path).strip():
        raise ValueError(f"Configure the {description} input path.")
    try:
        return pd.read_csv(path, dtype=str)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise ValueError(f"Cannot read the configured {description} CSV ({type(error).__name__}).") from None


def load_cohorts(args):
    """Load clinical metadata, expression features, and named partitions."""
    clinical = validate_cohort_table(
        _read_csv(args.clinical_csv, "clinical"),
        patient_id_column=getattr(args, "patient_id_column", None),
        source_labels=source_labels_from_args(args),
        participant_id_pattern=getattr(args, "participant_id_pattern", None),
        partition_labels=partition_labels_from_args(args),
    )
    rna = _read_csv(args.rna_csv, "RNA expression")
    if "sample_id" not in rna.columns:
        raise ValueError("RNA expression metadata requires a sample_id column.")
    if rna["sample_id"].isna().any():
        raise ValueError("RNA sample identifiers must not be missing.")
    rna["sample_id"] = rna["sample_id"].astype(str).str.strip()
    if (rna["sample_id"] == "").any() or rna["sample_id"].map(_identifier).duplicated().any():
        raise ValueError("RNA sample identifiers must be nonempty and globally unique.")
    feature_columns = [column for column in rna.columns if column != "sample_id"]
    if len(feature_columns) != int(args.RNAseq):
        raise ValueError("The RNA feature count does not match RNAseq in the configuration.")
    try:
        rna_values = rna[feature_columns].to_numpy(dtype=np.float32)
    except (TypeError, ValueError):
        raise ValueError("RNA features must be numeric.") from None
    if not np.isfinite(rna_values).all():
        raise ValueError("RNA features contain nonfinite values.")
    rna[feature_columns] = rna_values
    expected = set(clinical.loc[clinical["IS_SeqRNA"] == 1, "sample_id"])
    if not expected.issubset(set(rna["sample_id"])):
        raise ValueError("RNA features are missing for samples marked as RNA-available.")
    partitions = {
        "partition_01": ("SOURCE_A", "PARTITION_01"),
        "partition_02": ("SOURCE_A", "PARTITION_02"),
        "partition_03": ("SOURCE_A", "PARTITION_03"),
        "partition_04": ("SOURCE_B", "PARTITION_04"),
        "partition_05": ("SOURCE_C", "PARTITION_05"),
        "partition_06": ("SOURCE_D", "PARTITION_06"),
    }
    return {
        "clinical": clinical,
        "rna": rna,
        "partitions": {
            name: clinical.loc[(clinical["cohort"] == cohort) & (clinical["Strategy"] == strategy)].copy()
            for name, (cohort, strategy) in partitions.items()
        },
        "feature_directories": {
            cohort: getattr(args, f"wsi_feature_dir_{cohort.lower()}")
            for cohort in SOURCE_PARTITIONS
        },
        "partition_cohorts": {name: cohort for name, (cohort, _) in partitions.items()},
    }


def adjust_matrix(matrix, target_rows=3000):
    """Sample or zero-pad rows and mark padded indices with -1."""
    matrix = torch.as_tensor(matrix, dtype=torch.float32)
    if matrix.ndim != 2 or matrix.shape[0] == 0 or int(target_rows) <= 0:
        raise ValueError("WSI features require a nonempty matrix and positive patch count.")
    target_rows = int(target_rows)
    count = matrix.shape[0]
    if count > target_rows:
        indices = torch.randperm(count, device=matrix.device)[:target_rows]
        return matrix[indices], indices
    indices = torch.arange(count, device=matrix.device)
    if count < target_rows:
        matrix = torch.cat([matrix, matrix.new_zeros((target_rows - count, matrix.shape[1]))])
        indices = torch.cat([indices, indices.new_full((target_rows - count,), -1)])
    return matrix, indices


def _cache_settings(args):
    """Return settings that determine cached patch content and dimensions."""
    settings = {
        "version": _CACHE_VERSION,
        "max_patches": int(getattr(args, "wsi_cache_max_patches", 2000)),
        "seed": int(getattr(args, "seed", 42)),
        "feature_dim": int(args.D_WSI),
    }
    if settings["max_patches"] <= 0 or settings["feature_dim"] <= 0:
        raise ValueError("WSI cache patch count and feature dimension must be positive.")
    return settings


def _source_identity(datasets_path):
    """Hash a source location without persisting its private directory name."""
    source = os.path.normcase(str(Path(datasets_path).resolve()))
    return sha256(source.encode("utf-8")).hexdigest()


def _wsi_cache_path(datasets_path, args, sample_ids=None):
    """Derive a generic cache filename from source identity and settings."""
    settings = _cache_settings(args)
    selection = "all" if sample_ids is None else repr(sorted(set(sample_ids)))
    digest = sha256(f"{_source_identity(datasets_path)}|{sorted(settings.items())}|{selection}".encode("utf-8")).hexdigest()
    return Path(getattr(args, "wsi_cache_dir", "cache/wsi_features")) / f"wsi_{digest[:24]}.pt"


def _source_manifest(datasets_path, sample_ids=None):
    """Fingerprint source identity and HDF5 filenames, sizes, and timestamps."""
    try:
        directory = Path(datasets_path)
        if not directory.is_dir():
            raise ValueError("The configured WSI feature directory is unavailable.")
        files = sorted(directory.glob("*.h5"), key=lambda path: path.name.casefold())
        if sample_ids is not None:
            selected_ids = set(sample_ids)
            files = [path for path in files if path.stem in selected_ids]
            if {path.stem for path in files} != selected_ids:
                raise ValueError("WSI source files are missing for samples marked as WSI-available.")
        if not files:
            raise ValueError("The configured WSI feature directory contains no HDF5 files.")
        digest = sha256(_source_identity(directory).encode("ascii"))
        names = set()
        for path in files:
            name = _identifier(path.stem)
            if name in names:
                raise ValueError("WSI source files contain duplicate sample identifiers.")
            names.add(name)
            stat = path.stat()
            digest.update(f"{path.name}|{stat.st_size}|{stat.st_mtime_ns}\n".encode("utf-8"))
        return files, digest.hexdigest()
    except OSError as error:
        raise ValueError(f"Cannot inspect the configured WSI source ({type(error).__name__}).") from None


def _sample_patches(features, max_patches, rng):
    """Select a reproducible, order-preserving subset of existing patches."""
    if features.shape[0] <= max_patches:
        return features
    indices = np.sort(rng.choice(features.shape[0], size=max_patches, replace=False))
    return features[indices]


def _sample_rng(sample_id, seed):
    """Keep each sample's patch subset independent of other cohort records."""
    digest = sha256(f"{seed}|{sample_id}".encode("utf-8")).digest()
    return np.random.RandomState(int.from_bytes(digest[:4], "little"))


def _read_h5_features(path, feature_dim):
    """Read finite, nonempty embeddings from the HDF5 features dataset."""
    try:
        with h5py.File(path, "r") as handle:
            features = np.asarray(handle["features"][:], dtype=np.float32)
    except (OSError, KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Cannot read a required WSI features dataset ({type(error).__name__}).") from None
    if features.ndim != 2 or features.shape[0] == 0 or features.shape[1] != feature_dim:
        raise ValueError("WSI feature shape does not match the configured feature dimension.")
    if not np.isfinite(features).all():
        raise ValueError("WSI features contain nonfinite values.")
    return features


def build_wsi_feature_cache(datasets_path, args, cache_path=None, sample_ids=None):
    """Persist deterministic per-sample patches with path-free provenance.

    Caching fits no statistics and never combines samples. Fingerprints
    invalidate caches when locations, file metadata, or patch settings change.
    """
    settings = _cache_settings(args)
    files, fingerprint = _source_manifest(datasets_path, sample_ids)
    cache_path = Path(cache_path) if cache_path is not None else _wsi_cache_path(datasets_path, args, sample_ids)
    features = {}
    for path in files:
        values = _read_h5_features(path, settings["feature_dim"])
        values = _sample_patches(values, settings["max_patches"], _sample_rng(path.stem, settings["seed"]))
        # Float32 retains embedding precision and avoids float16 overflow.
        features[path.stem] = torch.from_numpy(np.ascontiguousarray(values)).clone()
    _, final_fingerprint = _source_manifest(datasets_path, sample_ids)
    if final_fingerprint != fingerprint:
        raise ValueError("WSI source files change during cache construction; retry with stable inputs.")
    payload = {"settings": settings, "source_fingerprint": fingerprint, "features": features}
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(dir=cache_path.parent, prefix="wsi_", suffix=".tmp", delete=False) as temporary:
            temporary_path = Path(temporary.name)
        torch.save(payload, temporary_path)
        temporary_path.replace(cache_path)
    except (OSError, RuntimeError) as error:
        raise ValueError(f"Cannot save the configured WSI cache ({type(error).__name__}).") from None
    print(f"[WSI cache] Saved {len(features)} samples; source location is masked.")
    return payload


def _cache_matches(payload, settings, fingerprint, files):
    """Validate cache provenance and tensor shapes before reusing embeddings."""
    if not isinstance(payload, dict) or payload.get("settings") != settings:
        return False
    if payload.get("source_fingerprint") != fingerprint:
        return False
    features = payload.get("features")
    if not isinstance(features, dict) or set(features) != {path.stem for path in files}:
        return False
    return all(
        torch.is_tensor(value) and value.ndim == 2
        and 0 < value.shape[0] <= settings["max_patches"]
        and value.shape[1] == settings["feature_dim"]
        and bool(torch.isfinite(value).all())
        for value in features.values()
    )


def get_or_build_wsi_cache(datasets_path, args, sample_ids=None):
    """Load a matching cache or rebuild it using tensor-only deserialization."""
    settings = _cache_settings(args)
    files, fingerprint = _source_manifest(datasets_path, sample_ids)
    key = (fingerprint, tuple(sorted(settings.items())))
    if key in _WSI_FEATURE_CACHE:
        return _WSI_FEATURE_CACHE[key]
    cache_path = _wsi_cache_path(datasets_path, args, sample_ids)
    payload = None
    if cache_path.exists() and not bool(getattr(args, "rebuild_wsi_cache", False)):
        try:
            payload = torch.load(cache_path, map_location="cpu", weights_only=True)
        except Exception:
            # Unreadable or incompatible caches rebuild from their source files.
            payload = None
    if not _cache_matches(payload, settings, fingerprint, files):
        payload = build_wsi_feature_cache(datasets_path, args, cache_path, sample_ids)
    _WSI_FEATURE_CACHE[key] = payload["features"]
    return payload["features"]


class MirrorDataset(Dataset):
    """Provide multimodal records with fixed-size patch sampling."""

    _sample_patches = True

    def __init__(self, clinical, rna_features, datasets_path, args):
        """Validate partition metadata and defer WSI reads until item access."""
        self.args = args
        self.datasets_path = datasets_path
        self.wsi_cache = None
        self._cache_initialized = False
        self.df = clinical.copy().reset_index(drop=True)
        if not self.df.empty:
            self.df = validate_cohort_table(
                self.df, patient_id_column=getattr(args, "patient_id_column", None),
                source_labels=source_labels_from_args(args),
                participant_id_pattern=getattr(args, "participant_id_pattern", None),
                partition_labels=partition_labels_from_args(args),
            )
        if self._sample_patches:
            if self.df.empty:
                raise ValueError("The requested optimization partition is empty.")
            if not ((self.df["cohort"] == "SOURCE_A") & (self.df["Strategy"] == "PARTITION_01")).all():
                raise ValueError("The optimization dataset contains an invalid partition.")
            if self.df["HRD_status"].isna().any():
                raise ValueError("Every optimization record requires a binary HRD label.")
        elif not self.df.empty and len(self.df[["cohort", "Strategy"]].drop_duplicates()) != 1:
            raise ValueError("Each dataset must contain one source and partition.")
        if int(args.N) <= 0 or int(args.D_WSI) <= 0 or int(args.RNAseq) <= 0:
            raise ValueError("Configured feature dimensions and patch counts must be positive.")
        if "sample_id" not in rna_features.columns:
            raise ValueError("RNA expression metadata requires a sample_id column.")
        if rna_features["sample_id"].map(_identifier).duplicated().any():
            raise ValueError("RNA sample identifiers must be unique.")
        # Each dataset retains RNA rows belonging to its own clinical partition.
        partition_ids = set(self.df["sample_id"])
        self.rna_features = rna_features.loc[
            rna_features["sample_id"].isin(partition_ids)
        ].set_index("sample_id")
        if self.rna_features.shape[1] != int(args.RNAseq):
            raise ValueError("The RNA feature count does not match RNAseq in the configuration.")

    def __len__(self):
        """Return the number of records in this partition."""
        return len(self.df)

    def _wsi_features(self, sample_id):
        """Read a sample lazily and supply a mask for genuine patch rows."""
        settings = _cache_settings(self.args)
        if bool(getattr(self.args, "use_wsi_cache", True)):
            if not self._cache_initialized:
                sample_ids = self.df.loc[self.df["IS_FFPE"] == 1, "sample_id"].tolist()
                self.wsi_cache = get_or_build_wsi_cache(self.datasets_path, self.args, sample_ids)
                self._cache_initialized = True
            if sample_id not in self.wsi_cache:
                raise ValueError("WSI features are missing for a sample marked as WSI-available.")
            features = self.wsi_cache[sample_id]
        else:
            features = _read_h5_features(Path(self.datasets_path) / f"{sample_id}.h5", settings["feature_dim"])
            features = _sample_patches(features, settings["max_patches"], _sample_rng(sample_id, settings["seed"]))
        if self._sample_patches:
            features, indices = adjust_matrix(features, self.args.N)
            return features, indices >= 0
        features = torch.as_tensor(features, dtype=torch.float32)
        return features, torch.ones(features.shape[0], dtype=torch.bool)

    def __getitem__(self, index):
        """Return features, modality availability, patch masks, and outcomes."""
        row = self.df.iloc[index]
        sample_id = row["sample_id"]
        has_wsi = bool(row["IS_FFPE"])
        has_rna = bool(row["IS_SeqRNA"])
        if has_wsi:
            wsi, pad_mask = self._wsi_features(sample_id)
        else:
            wsi = torch.zeros(int(self.args.N), int(self.args.D_WSI), dtype=torch.float32)
            pad_mask = torch.zeros(int(self.args.N), dtype=torch.bool)
        if has_rna:
            if sample_id not in self.rna_features.index:
                raise ValueError("RNA features are missing for a sample marked as RNA-available.")
            try:
                rna = torch.as_tensor(self.rna_features.loc[sample_id].to_numpy(dtype=np.float32))
            except (TypeError, ValueError):
                raise ValueError("RNA features must be numeric.") from None
            if not bool(torch.isfinite(rna).all()):
                raise ValueError("RNA features contain nonfinite values.")
        else:
            rna = torch.zeros(int(self.args.RNAseq), dtype=torch.float32)
        return {
            "index": index,
            "WSI_feature": wsi,
            "RNASeq_feature": rna,
            "pad_mask": pad_mask,
            "mask_wsi": torch.tensor(has_wsi, dtype=torch.long),
            "mask_SeqRNA": torch.tensor(has_rna, dtype=torch.long),
            "patient_name": sample_id,
            "OS": torch.tensor(row["OStime"], dtype=torch.float32),
            "OSState": torch.tensor(row["OS"], dtype=torch.float32),
            "Label": torch.tensor(row["HRD_status"], dtype=torch.float32),
            "DFS": torch.tensor(row["DFStime"], dtype=torch.float32),
            "DFSState": torch.tensor(row["DFS"], dtype=torch.float32),
        }


class MirrorRecordDataset(MirrorDataset):
    """Provide deterministic, bounded patch bags for individual records."""

    _sample_patches = False


def collate_multimodal(batch):
    """Stack fixed-shape tensors and retain sample identifiers."""
    if not batch:
        raise ValueError("Cannot collate an empty multimodal batch.")
    result = {}
    for key in batch[0]:
        values = [sample[key] for sample in batch]
        result[key] = torch.stack(values) if torch.is_tensor(values[0]) else values
    return result
