"""Generate patient-grouped clinical partition assignments."""

import argparse
from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split

from .multimodal_dataset import (
    _read_csv,
    patient_group_keys,
    source_labels_from_args,
    validate_cohort_table,
)


def split_cohorts(frame, partition_fraction=0.7, seed=2026, patient_id_column=None,
                  source_labels=None, participant_id_pattern=None):
    """Assign each patient's specimens together with patient-level stratification."""
    if not 0 < partition_fraction < 1:
        raise ValueError("The partition fraction must lie strictly between zero and one.")
    frame = validate_cohort_table(
        frame, patient_id_column, require_strategy=False,
        source_labels=source_labels, participant_id_pattern=participant_id_pattern,
    )
    frame["Strategy"] = pd.NA
    for source, partition in (("SOURCE_B", "PARTITION_04"),
                              ("SOURCE_C", "PARTITION_05"),
                              ("SOURCE_D", "PARTITION_06")):
        frame.loc[frame["cohort"] == source, "Strategy"] = partition
    source_a = frame.loc[frame["cohort"] == "SOURCE_A"].copy()
    if source_a.empty:
        raise ValueError("Split generation requires SOURCE_A records.")
    source_a["_patient_group"] = patient_group_keys(
        source_a, patient_id_column, participant_id_pattern
    )
    grouped = source_a.groupby("_patient_group", sort=True)
    if (grouped["HRD_status"].nunique(dropna=True) > 1).any():
        raise ValueError("A SOURCE_A patient has conflicting HRD labels across specimens.")
    labeled_count = grouped["HRD_status"].count()
    if ((labeled_count > 0) & (labeled_count < grouped.size())).any():
        raise ValueError("A SOURCE_A patient mixes labeled and unlabeled specimens; reconcile metadata before splitting.")
    patients = grouped.agg(HRD_status=("HRD_status", "first"), IS_FFPE=("IS_FFPE", "max"))
    labeled_patients = patients.loc[patients["HRD_status"].notna()]
    unlabeled_groups = set(patients.index[patients["HRD_status"].isna()])
    frame.loc[source_a.index[source_a["_patient_group"].isin(unlabeled_groups)], "Strategy"] = "PARTITION_03"
    if len(labeled_patients) < 4 or labeled_patients["HRD_status"].nunique() != 2:
        raise ValueError("Stratification requires enough observations from both HRD classes.")

    # Patient-level stratification keeps related specimens in one partition.
    combined = labeled_patients["IS_FFPE"].astype(int).astype(str) + "_" + labeled_patients["HRD_status"].astype(int).astype(str)
    n_first = int(partition_fraction * len(labeled_patients))
    n_second = len(labeled_patients) - n_first
    strata = combined
    if combined.value_counts().min() < 2 or combined.nunique() > min(n_first, n_second):
        strata = labeled_patients["HRD_status"]
    try:
        first_groups, second_groups = train_test_split(
            labeled_patients.index, train_size=partition_fraction,
            stratify=strata, random_state=seed, shuffle=True,
        )
    except ValueError:
        raise ValueError("Patient counts do not support the requested stratified partition sizes.") from None
    frame.loc[source_a.index[source_a["_patient_group"].isin(first_groups)], "Strategy"] = "PARTITION_01"
    frame.loc[source_a.index[source_a["_patient_group"].isin(second_groups)], "Strategy"] = "PARTITION_02"
    return validate_cohort_table(
        frame, patient_id_column, participant_id_pattern=participant_id_pattern
    )


def summarize_partitions(frame):
    """Summarize labels and available modalities for each assigned partition."""
    rows = []
    order = ("PARTITION_01", "PARTITION_02", "PARTITION_03", "PARTITION_04", "PARTITION_05", "PARTITION_06")
    for strategy in order:
        subset = frame.loc[frame["Strategy"] == strategy]
        if subset.empty:
            continue
        has_rna = subset["IS_SeqRNA"] == 1
        has_wsi = subset["IS_FFPE"] == 1
        rows.append({
            "Strategy": strategy,
            "samples": len(subset),
            "HRD_positive": int((subset["HRD_status"] == 1).sum()),
            "HRD_negative": int((subset["HRD_status"] == 0).sum()),
            "paired_modalities": int((has_rna & has_wsi).sum()),
            "RNA_only": int((has_rna & ~has_wsi).sum()),
            "WSI_only": int((has_wsi & ~has_rna).sum()),
            "OS_observed": int((subset["OS"].notna() & subset["OStime"].notna()).sum()),
        })
    return pd.DataFrame(rows)


def main():
    """Generate explicitly requested split and summary CSV files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input clinical CSV.")
    parser.add_argument("--output", required=True, type=Path, help="Clinical CSV with partition assignments.")
    parser.add_argument("--summary", type=Path, help="Optional partition summary CSV.")
    parser.add_argument("--config", type=Path, help="Optional local YAML configuration.")
    parser.add_argument("--partition-fraction", type=float, default=0.7)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--patient-id-column", default=None)
    parser.add_argument("--participant-id-pattern", default=None)
    args = parser.parse_args()
    from utils.config import parse_config

    settings = parse_config(["--config", str(args.config)] if args.config else [])
    patient_column = args.patient_id_column or getattr(settings, "patient_id_column", None)
    participant_pattern = args.participant_id_pattern or getattr(settings, "participant_id_pattern", None)
    frame = split_cohorts(
        _read_csv(args.input, "clinical"), args.partition_fraction, args.seed,
        patient_column, source_labels=source_labels_from_args(settings),
        participant_id_pattern=participant_pattern,
    )
    summary = summarize_partitions(frame)
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, index=False)
        if args.summary:
            args.summary.parent.mkdir(parents=True, exist_ok=True)
            summary.to_csv(args.summary, index=False)
    except OSError as error:
        raise ValueError(f"Cannot write the requested split outputs ({type(error).__name__}); paths are masked.") from None
    print(summary.to_string(index=False))
    print("Partition assignments are saved to the configured output location.")


if __name__ == "__main__":
    main()
