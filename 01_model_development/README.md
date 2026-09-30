# MIRROR

MIRROR is a multimodal framework for predicting homologous recombination deficiency (HRD) from whole-slide image (WSI) embeddings and RNA expression features. It combines gated attention pooling, an RNA encoder, bidirectional cross-attention, and learned representations for unavailable modalities.

The repository covers three major stages:

1. **Data preparation:** clinical metadata checks, patient grouping, HDF5 inspection, and deterministic patch-feature caching.
2. **Model development:** multimodal feature encoding, objective computation, parameter optimization, and checkpoint management.
3. **Output analysis:** HRD prediction, classification and survival metrics, and structured result export.

The codebase provides a configurable research workflow with separate modules for data handling, neural components, numerical analysis, and persistence.

## Repository Layout

```text
MIRROR/
|-- config/
|   `-- config.yaml
|-- datasets/
|   |-- __init__.py
|   |-- data_loaders.py
|   |-- inspect_hdf5.py
|   |-- multimodal_dataset.py
|   `-- split_cohorts.py
|-- models/
|   |-- __init__.py
|   |-- attention.py
|   |-- augmentation.py
|   |-- contrastive_loss.py
|   |-- mirror.py
|   |-- resnet.py
|   `-- survival.py
|-- utils/
|   |-- __init__.py
|   |-- config.py
|   |-- metric_log.py
|   |-- save_model.py
|   |-- survival.py
|   `-- yaml_config_hook.py
|-- build_wsi_cache.py
|-- evaluation.py
|-- optimize.py
|-- requirements.txt
|-- .gitignore
`-- README.md
```

## Important Note

All paths and source identifiers in the public configuration are generic placeholders. Actual storage locations and metadata-label values belong in a local configuration file, such as `config/experiment.local.yaml`. Local overrides and the default input, cache, and output directories are excluded from source control. Add matching ignore rules for custom storage locations.

This repository consumes pre-extracted WSI embeddings and prepared RNA expression vectors. Raw-slide tiling, foundation-model feature extraction, and expression preprocessing are not included.

## Environment

The pipeline requires a Python environment with the packages listed in `requirements.txt`. Core dependencies include:

- `torch`
- `numpy`
- `pandas`
- `h5py`
- `PyYAML`
- `scikit-learn`
- `lifelines`
- `python-docx`

Optional image and encoder utilities also use `torchvision`, `opencv-python`, and `Pillow`.

```bash
pip install -r requirements.txt
```

Package versions are not pinned. Compatibility with a particular Python, PyTorch, or CUDA environment requires separate verification.

## Workflow

### 01. Data Preparation

#### 01.1 Organize input files

Prepare a clinical metadata table, an RNA expression matrix, and directories of slide-level feature files.

```text
data/
|-- clinical.csv
|-- rna_expression.csv
`-- wsi/
    |-- source_a/
    |   |-- sample_001.h5
    |   `-- sample_002.h5
    |-- source_b/
    |-- source_c/
    `-- source_d/
```

The clinical table contains the following columns:

```text
sample_id,cohort,Strategy,IS_SeqRNA,IS_FFPE,HRD_status,OS,OStime,DFS,DFStime
sample_001,SOURCE_A,PARTITION_01,1,1,1,0,36,0,30
sample_002,SOURCE_A,PARTITION_02,1,1,0,1,24,1,18
sample_003,SOURCE_B,PARTITION_04,0,1,1,0,42,0,35
```

These rows are illustrative. Availability flags use 0 or 1, and each observation requires at least one available modality. Binary outcomes accept 0, 1, or missing values. Recorded durations use a consistent unit.

The RNA matrix stores one sample per row and one feature per column:

```text
sample_id,gene_0001,gene_0002,gene_0003,...
sample_001,0.1532,-0.2910,0.0084,...
sample_002,0.1321,-0.1442,0.0213,...
```

The number of feature columns matches `RNAseq` in the configuration. Feature order remains consistent across all inputs.

#### 01.2 Inspect WSI feature files

Each `<sample_id>.h5` file contains a `features` dataset with shape `(patch_count, D_WSI)`. Feature values are finite numeric embeddings.

```bash
python -m datasets.inspect_hdf5 data/wsi/source_a/sample_001.h5
```

The inspection utility reports object names, array shapes, and data types without loading complete feature arrays.

#### 01.3 Configure metadata and storage

Use `config/config.yaml` for public defaults and an ignored local YAML file for machine-specific values. A local override can contain a subset of settings:

```yaml
clinical_csv: data/clinical.csv
rna_csv: data/rna_expression.csv
wsi_feature_dir_source_a: data/wsi/source_a
wsi_feature_dir_source_b: data/wsi/source_b
wsi_feature_dir_source_c: data/wsi/source_c
wsi_feature_dir_source_d: data/wsi/source_d
out_dir: results

D_WSI: 1024
RNAseq: 12042
N: 500
batch_size: 4
```

The `source_*_label` and `partition_*_label` settings map metadata values to neutral internal identifiers. Relative paths in this configuration resolve from the project root.

An explicit `patient_id_column` or a `participant_id_pattern` with a capture group supports matching related specimens. Identical sample names and conflicting patient assignments produce an error.

#### 01.4 Optionally generate partition assignments

The metadata utility creates patient-grouped assignments and an optional summary table.

```bash
python -m datasets.split_cohorts --input data/clinical_unsplit.csv --output data/clinical.csv --summary data/partition_summary.csv --config config/experiment.local.yaml
```

Use a separate output file to retain the original metadata. This utility is an explicit command and does not run automatically when the main entry point starts.

#### 01.5 Prepare WSI caches

The cache utility stores deterministic patch subsets and checks source fingerprints before reusing existing files.

```bash
python build_wsi_cache.py --config config/experiment.local.yaml
```

To refresh cached feature files:

```bash
python build_wsi_cache.py --config config/experiment.local.yaml --rebuild_wsi_cache true
```

Cache filenames and provenance use digests rather than plaintext source directories. The cache stores float32 embeddings; `wsi_cache_max_patches` controls the maximum retained patch count.

### 02. Model Development

#### 02.1 Configure MIRROR

The model uses a gated attention module to summarize WSI patch embeddings and a multilayer RNA encoder to represent expression vectors. Bidirectional cross-attention combines the two representations. Learned tokens represent unavailable modalities, while patch masks distinguish observed embeddings from padding.

The main settings include:

- `D_WSI`, `RNAseq`, and `N`: feature dimensions and patch-bag size.
- `start_lr`, `weight_decay`, and `Epoch`: optimizer and scheduling settings.
- `w_cls`, `w_uni`, and `w_align`: classification, auxiliary, and alignment weights.
- `eval_every`: reporting and checkpoint interval.
- `use_wsi_cache` and `wsi_cache_max_patches`: feature-cache settings.

#### 02.2 Run the optimization workflow

```bash
python optimize.py --config config/experiment.local.yaml
```

Command-line options override YAML values. Boolean options accept explicit values:

```bash
python optimize.py --config config/experiment.local.yaml --use_wsi_cache false
```

Each run writes to a separate timestamped directory under the configured output location.

#### 02.3 Continue from a saved checkpoint

Set `initial_checkpoint` to a resumable checkpoint in the local configuration:

```yaml
initial_checkpoint: results/MIRROR/run_identifier/last_model.pth
```

Continuation restores parameter, optimizer, scheduler, epoch, and random-generator state. Version markers, sample fingerprints, and the configured epoch horizon must match. Keep feature content and gene order unchanged.

### 03. Output Analysis

#### 03.1 Inspect prediction tables

The workflow exports `predictions.csv` with the following fields:

```text
sample_id,HRD_score,HRD_prob,PRE_HRD,Label,OStime,OS,dataset,modality
```

`HRD_score` is the model logit. `HRD_prob` is its sigmoid probability, and `PRE_HRD` is the binary assignment at the corresponding threshold. The `modality` field distinguishes Multimodal, WSI, and RNASeq outputs. The `dataset` field uses neutral identifiers such as `set_01`.

#### 03.2 Review numerical summaries

`evaluation_metrics.csv` contains classification accuracy, AUROC, sensitivity, specificity, concordance, log-rank significance, and a Cox hazard ratio where computable. A value of -1 denotes an unavailable metric.

The concordance calculation uses the convention that a higher HRD score corresponds to longer predicted survival. Survival statistics describe associations in the supplied observations.

#### 03.3 Locate saved artifacts

```text
results/
`-- MIRROR/
    `-- run_identifier/
        |-- best_model.pth
        |-- last_model.pth
        |-- mirror_model.pth
        |-- thresholds.json
        |-- predictions.csv
        |-- evaluation_metrics.csv
        `-- run_metrics.docx
```

`mirror_model.pth` contains the selected parameters, model dimensions, and thresholds for downstream use. It does not contain optimizer continuation state. Numerical thresholds are also available in `thresholds.json`.

`best_model.pth` is written when the comparator accepts a new candidate; `last_model.pth` records the most recent scheduled state.

## Module Guide

- **`datasets/multimodal_dataset.py`**  
  Checks metadata, resolves neutral identifiers, loads multimodal features, and manages deterministic caches.
- **`datasets/data_loaders.py`**  
  Constructs a batch loader and named dataset objects with reproducible sampling.
- **`datasets/split_cohorts.py`**  
  Generates patient-grouped assignments and summarizes labels and modality availability.
- **`models/mirror.py`**  
  Implements MIRROR encoders, missing-modality representations, cross-attention fusion, and objective components.
- **`optimize.py`**  
  Coordinates parameter updates, checkpoint management, state restoration, and output generation.
- **`evaluation.py`**  
  Collects scores, calculates thresholds and metrics, and exports structured tables.
- **`utils/config.py`**  
  Combines public defaults, local overrides, and command-line options, with masking for exported private settings.
- **`utils/metric_log.py`**  
  Records messages and numerical summaries in the console and a Word document.

## Review Status

The source receives static inspection only. Python and project entry points are not executed during this code-editing review. Runtime behavior and dependency compatibility remain unverified.

## Citation

Citation details can be added when a publication DOI or preprint identifier is available.
