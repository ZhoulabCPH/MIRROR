"""Prepare and validate WSI feature caches."""

from utils.config import parse_config


def main(argv=None):
    """Build caches for the configured sample groups."""
    args = parse_config(argv, description=__doc__)
    from datasets.multimodal_dataset import get_or_build_wsi_cache, load_cohorts

    cohorts = load_cohorts(args)
    for index, (partition, clinical) in enumerate(cohorts["partitions"].items(), start=1):
        sample_ids = clinical.loc[clinical["IS_FFPE"] == 1, "sample_id"].tolist()
        if not sample_ids:
            continue
        cohort = cohorts["partition_cohorts"][partition]
        feature_dir = cohorts["feature_directories"][cohort]
        print(f"[WSI cache] Preparing set_{index:02d}: {len(sample_ids)} slides.")
        get_or_build_wsi_cache(feature_dir, args, sample_ids=sample_ids)


if __name__ == "__main__":
    main()
