"""Inspect an explicitly supplied HDF5 file without loading feature arrays."""

import argparse
from pathlib import Path

import h5py


def inspect_hdf5(file_path):
    """Print dataset names, shapes, and dtypes while masking the source path."""
    try:
        with h5py.File(Path(file_path), mode="r") as handle:
            print("Root objects:", list(handle.keys()))

            def show_structure(name, item):
                """Display structural metadata for each visited HDF5 object."""
                if isinstance(item, h5py.Dataset):
                    print(f"Dataset: {name}, shape={item.shape}, dtype={item.dtype}")
                elif isinstance(item, h5py.Group):
                    print(f"Group: {name}")

            handle.visititems(show_structure)
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot inspect the supplied HDF5 file ({type(error).__name__}); its path is masked.") from None


def main():
    """Parse the file argument and report structure without import-time I/O."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="Local HDF5 feature file to inspect.")
    args = parser.parse_args()
    inspect_hdf5(args.file)


if __name__ == "__main__":
    main()
