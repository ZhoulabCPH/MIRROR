"""Parse MIRROR settings and redact private values from exports."""

import argparse
import math
import os
from pathlib import Path

from .yaml_config_hook import yaml_config_hook


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "config.yaml"
SOURCE_KEYS = ("source_a", "source_b", "source_c", "source_d")
PARTITION_KEYS = tuple(f"partition_{index:02d}" for index in range(1, 7))
PATH_KEYS = frozenset({
    "clinical_csv", "rna_csv", "out_dir", "wsi_cache_dir",
    "initial_checkpoint", "wsi_dir", "config",
    *(f"wsi_feature_dir_{source}" for source in SOURCE_KEYS),
    *(f"wsi_path_{source}" for source in SOURCE_KEYS),
})
OPTIONAL_PATH_KEYS = frozenset({
    "initial_checkpoint", "wsi_dir",
    *(f"wsi_path_{source}" for source in SOURCE_KEYS),
})
PRIVATE_KEYS = PATH_KEYS | frozenset({
    "patient_id_column", "participant_id_pattern",
    *(f"{source}_label" for source in SOURCE_KEYS),
    *(f"{partition}_label" for partition in PARTITION_KEYS),
})
INTEGER_KEYS = frozenset({
    "seed", "Epoch", "batch_size", "NUM_CLASSES", "D_WSI", "N", "RNAseq",
    "workers", "eval_every", "wsi_cache_max_patches",
})
FLOAT_KEYS = frozenset({"start_lr", "weight_decay", "w_cls", "w_uni", "w_align"})
BOOLEAN_KEYS = frozenset({"use_wsi_cache", "rebuild_wsi_cache"})


def parse_bool(value):
    """Interpret explicit Boolean values."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "on"}:
        return True
    if normalized in {"false", "0", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError("Use true or false for a Boolean option.")


def _value_type(key, value):
    """Choose the scalar type for a configuration field."""
    if key in BOOLEAN_KEYS or isinstance(value, bool):
        return parse_bool
    if key in INTEGER_KEYS:
        return int
    if key in FLOAT_KEYS:
        return float
    return str if value is None else type(value)


def _resolve_path(value):
    """Resolve a configured path or retain an unset optional value."""
    if value is None or str(value).strip().lower() in {"", "none", "null"}:
        return None
    path = Path(os.path.expandvars(str(value))).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return str(path.resolve())


def mask_path(value):
    """Return a generic display value for a storage location."""
    return None if value is None else "<path>"


def public_config(args):
    """Return settings with storage and source-identification fields masked."""
    settings = vars(args) if isinstance(args, argparse.Namespace) else dict(args)
    return {
        key: (None if value is None else f"<{key}>") if key in PRIVATE_KEYS else value
        for key, value in settings.items()
    }


def parse_config(argv=None, description=None):
    """Merge public defaults, a local YAML override, and typed CLI options."""
    preliminary = argparse.ArgumentParser(add_help=False)
    preliminary.add_argument("--config", default=str(DEFAULT_CONFIG))
    selected, _ = preliminary.parse_known_args(argv)
    selected_path = Path(selected.config).expanduser()
    if not selected_path.is_absolute():
        selected_path = PROJECT_ROOT / selected_path

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", default=str(selected_path), help="YAML override file.")
    try:
        config = yaml_config_hook(DEFAULT_CONFIG)
        if selected_path.resolve() != DEFAULT_CONFIG.resolve():
            overrides = yaml_config_hook(selected_path)
            if set(overrides).difference(config):
                parser.error("The override file contains an unsupported configuration key.")
            config.update(overrides)
    except (OSError, ValueError) as error:
        detail = str(error) if not isinstance(error, OSError) else "The selected file cannot be read."
        parser.error(f"Cannot load configuration. {detail}")

    for key, value in config.items():
        if not isinstance(key, str) or not key.replace("_", "").isalnum():
            parser.error("Configuration keys must contain letters, digits, and underscores.")
        if key == "config":
            parser.error("The config key is reserved for the command-line file selector.")
        if isinstance(value, (list, dict)):
            parser.error(f"Configuration option {key} must contain a scalar value.")
        if key in INTEGER_KEYS | FLOAT_KEYS | BOOLEAN_KEYS and value is None:
            parser.error(f"Configuration option {key} cannot be null.")
        if key in INTEGER_KEYS and isinstance(value, float) and not value.is_integer():
            parser.error(f"Configuration option {key} must contain an integer.")
        option_type = _value_type(key, value)
        try:
            default = None if value is None else option_type(value)
        except (ValueError, TypeError, argparse.ArgumentTypeError):
            parser.error(f"Invalid value for configuration option {key}.")
        flags = [f"--{key}"]
        if "_" in key:
            flags.append(f"--{key.replace('_', '-')}")
        if option_type is parse_bool:
            parser.add_argument(*flags, dest=key, default=default, nargs="?", const=True, type=parse_bool)
            parser.add_argument(f"--no-{key.replace('_', '-')}", dest=key, action="store_false")
        else:
            parser.add_argument(*flags, dest=key, default=default, type=option_type)

    args = parser.parse_args(argv)
    for key in PATH_KEYS:
        if hasattr(args, key):
            setattr(args, key, _resolve_path(getattr(args, key)))
    for key in INTEGER_KEYS - {"seed", "workers"}:
        if getattr(args, key, 1) <= 0:
            parser.error(f"{key} must be positive.")
    if args.workers < 0 or args.seed < 0:
        parser.error("workers and seed must be nonnegative.")
    if any(not math.isfinite(getattr(args, key)) for key in FLOAT_KEYS):
        parser.error("Learning rate, weight decay, and loss weights must be finite.")
    if args.start_lr <= 0 or any(getattr(args, key) < 0 for key in FLOAT_KEYS - {"start_lr"}):
        parser.error("Learning rate must be positive; loss weights and weight decay must be nonnegative.")
    for key in PATH_KEYS - OPTIONAL_PATH_KEYS:
        if hasattr(args, key) and getattr(args, key) is None:
            parser.error(f"{key} must specify a path.")
    for keys in (SOURCE_KEYS, PARTITION_KEYS):
        labels = [getattr(args, f"{key}_label") for key in keys]
        if any(label is None or not str(label).strip() for label in labels):
            parser.error("Labels must be nonempty.")
        if len({str(label).strip().casefold() for label in labels}) != len(labels):
            parser.error("Labels within each group must be distinct.")
    return args
