"""Load safe YAML mappings and merge optional configuration defaults."""

from pathlib import Path

import yaml


def yaml_config_hook(config_file, _active_files=None):
    """Merge included defaults, then apply values from the requested file.

    Each ``defaults`` entry contains one directory-to-filename mapping. Includes
    resolve relative to the containing YAML file and accept a ``.yaml`` suffix.
    The loader rejects recursive includes and never mutates parsed defaults.
    """
    config_path = Path(config_file).expanduser().resolve()
    active_files = set() if _active_files is None else set(_active_files)
    if config_path in active_files:
        raise ValueError("Configuration defaults contain a circular include.")
    active_files.add(config_path)
    try:
        with config_path.open("r", encoding="utf-8") as stream:
            config = yaml.safe_load(stream)
    except yaml.YAMLError:
        raise ValueError("The configuration contains invalid YAML.") from None
    if config is None:
        return {}
    if not isinstance(config, dict):
        raise ValueError("The configuration must contain a YAML mapping.")

    defaults = config.get("defaults", [])
    if not isinstance(defaults, list):
        raise ValueError("Configuration defaults must contain a list of mappings.")
    merged = {}
    for entry in defaults:
        if not isinstance(entry, dict) or len(entry) != 1:
            raise ValueError("Each configuration default must contain one mapping.")
        directory, filename = next(iter(entry.items()))
        if not isinstance(directory, str) or not isinstance(filename, str):
            raise ValueError("Configuration default paths must contain strings.")
        filename = filename if filename.endswith((".yaml", ".yml")) else filename + ".yaml"
        included_path = config_path.parent / directory / filename
        merged.update(yaml_config_hook(included_path, active_files))
    merged.update({key: value for key, value in config.items() if key != "defaults"})
    return merged
