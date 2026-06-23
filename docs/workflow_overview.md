````markdown
# TRACE-Omics workflow overview

**TRACE-Omics** stands for **Transcriptomic Robustness Assessment and Cross-cohort Evaluation for Omics machine learning**.

The repository implements a staged workflow for evaluating robust and generalisable transcriptomic machine learning using multi-cohort RNA-seq data. The central assumption is that reliable omics learning cannot be established from high internal performance alone. Instead, each stage evaluates whether disease-associated signal is preserved, stabilised, learned, transported, and biologically interpretable.

## Conceptual framework

TRACE-Omics evaluates transcriptomic learning as a sequence of evidence-generating stages:

```text
Representation quality
        ↓
Feature-space stability and conditioning
        ↓
Predictive benchmarking
        ↓
Nested cross-validation and optimisation stability
        ↓
Cross-cohort validation
        ↓
Statistical and biological reliability
````

Each stage asks a different question.

| Stage                      | Core question                                                                               | Main risk being tested                                            |
| -------------------------- | ------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Harmonisation diagnostics  | Does the representation reduce technical structure while preserving outcome signal?         | Batch effects, overcorrection, residual depth dependence          |
| Feature-space conditioning | Are selected genes stable, non-redundant and informative?                                   | Feature instability, redundancy, high-dimensional noise           |
| Benchmarking               | Can learners recover predictive signal across HC-FSF-MLA combinations?                      | Algorithm-specific overfitting or weak learnability               |
| Nested CV                  | Does performance remain stable after hyperparameter optimisation?                           | Optimisation instability and favourable split selection           |
| Cross-cohort validation    | Does performance transport across held-out cohorts?                                         | Distribution shift and cohort-specific shortcut learning          |
| Reliability analysis       | Do final models remain statistically and biologically interpretable on a held-out test set? | Poor calibration, unstable biomarkers and weak biological support |

## Repository structure

```text
trace-omics/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── config/
│   └── example_config.yml
├── data/
│   ├── README.md
│   ├── raw/
│   └── processed/
├── docs/
│   └── workflow_overview.md
├── R/
│   └── utils_io.R
├── scripts/
│   ├── 01_hc_representation_diagnostics.R
│   ├── 02_fsf_stability_conditioning.R
│   ├── 03_benchmarking_ncv_ccv.R
│   └── 04_statistical_biological_reliability.R
└── results/
```

Raw data, processed data and generated outputs are intentionally excluded from version control. The repository stores code, configuration files, documentation and reusable helper functions.

## Input data requirements

TRACE-Omics expects disease-specific data folders.

Example:

```text
data/raw/TB/
├── counts/
│   ├── cohort_1.csv
│   ├── cohort_2.csv
│   └── cohort_3.csv
└── metadata.rds
```

Each count table must contain:

| Requirement       | Description                                               |
| ----------------- | --------------------------------------------------------- |
| `symbol` column   | Gene symbol column                                        |
| Sample columns    | One column per RNA-seq sample                             |
| Count-like values | Raw or count-like expression values before transformation |

The metadata file must contain:

| Column      | Description                                         |
| ----------- | --------------------------------------------------- |
| `sample_id` | Sample identifier matching count-table column names |
| `cohort`    | Cohort, study, batch or dataset label               |
| `outcome`   | Binary class label, usually `case` or `control`     |

The scripts assume that `case` is the positive class and `control` is the negative class unless changed in the configuration or script header.

## Script 01: Harmonisation-configuration comparisons and representation diagnostics

**File:** `scripts/01_hc_representation_diagnostics.R`

This script evaluates whether the transcriptomic representation is suitable for downstream learning.

It compares:

```text
RAW counts
VST only
VST + ComBat
VST + ComBat-Ref
VST + SVA
ComBat-Seq + log1p
```

Main outputs:

```text
results/01_hc_representation_diagnostics/<disease>/
├── <disease>_harmonised_matrices.rds
├── <disease>_qc_summary.csv
├── <disease>_depth_dependence_overall.csv
├── <disease>_depth_dependence_by_cohort.csv
├── <disease>_pca_scores.csv
├── <disease>_pca_variance.csv
├── <disease>_silhouette_summary.csv
├── <disease>_qc_total_signal.png
├── <disease>_pca_hc_comparison.png
├── <disease>_silhouette_profiles.png
└── <disease>_sessionInfo.txt
```

Interpretation:

* Good harmonisation reduces cohort structure without removing disease-associated structure.
* Residual depth dependence suggests that transformed signal may still reflect sequencing depth.
* PCA and silhouette diagnostics separate visual structure from quantitative evidence.
* The selected representation should support downstream feature selection and modelling.

## Script 02: Feature-selection framework comparisons and feature-space conditioning

**File:** `scripts/02_fsf_stability_conditioning.R`

This script evaluates whether selected feature spaces are stable, informative and well-conditioned.

It compares feature-selection frameworks such as:

```text
ANOVA
AUC
RF_Imp
RF_Var
Boruta
CMIM
mRMR
Info_Gain
Hybrid A
Hybrid B
```

Main outputs:

```text
results/02_fsf_stability_conditioning/<disease>/
├── feature_subsets/
├── figures/
│   ├── <disease>_single_filter_JI_heatmap.png
│   ├── <disease>_single_filter_eRank_vs_meanVIF.png
│   └── <disease>_hybrid_filter_eRank_vs_meanVIF.png
└── tables/
    ├── <disease>_single_filter_jaccard_index.csv
    ├── <disease>_single_filter_stability_conditioning.csv
    ├── <disease>_hybrid_filter_stability_conditioning.csv
    └── <disease>_sessionInfo.txt
```

Interpretation:

* Jaccard Index evaluates feature-set overlap.
* Effective rank evaluates the intrinsic dimensionality of the selected feature space.
* VIF evaluates redundancy and collinearity.
* A useful feature space should balance stability, low redundancy and sufficient informational diversity.

## Script 03: Benchmarking, nested cross-validation and cross-cohort validation

**File:** `scripts/03_benchmarking_ncv_ccv.R`

This script evaluates predictive learnability, optimisation stability and transportability.

It includes:

1. Baseline benchmarking across HCs, FSFs and MLAs.
2. Prediction-level output extraction.
3. ROC and calibration plots for selected configurations.
4. Nested cross-validation with Bayesian tuning.
5. Hyperparameter sensitivity analysis.
6. Outer CV using best nested-CV parameters.
7. Cross-cohort validation.
8. CCV ROC, PR, calibration and permutation testing.

Main outputs:

```text
results/03_benchmarking_ncv_ccv/<disease>/
├── 01_baseline_benchmark/
├── 02_nested_cv/
├── 03_outercv_from_bestparams/
├── 04_cross_cohort_validation/
├── figures/
└── summary/
```

Interpretation:

* Baseline benchmarking identifies which learners can exploit the representation.
* Nested CV tests whether tuned performance is stable under internal validation.
* Hyperparameter sensitivity distinguishes broad stable regions from narrow optimisation peaks.
* CCV tests whether the model transports across cohort-level distribution shift.
* ROC, PR and calibration curves should be interpreted together, not as interchangeable summaries.

## Script 04: Statistical and biological reliability

**File:** `scripts/04_statistical_biological_reliability.R`

This script evaluates final selected models using a held-out test set drawn from the same expected distribution as the training data.

It includes:

1. Final model training.
2. Production bundle creation.
3. Held-out test prediction.
4. Youden threshold optimisation.
5. ROC, PR and calibration analysis.
6. DALEX permutation feature importance.
7. Consensus biomarker ranking.
8. Reactome pathway enrichment.

Main outputs:

```text
results/04_statistical_biological_reliability/<disease>/
├── production_bundles/
├── heldout_test/
├── feature_importance/
├── reactome_enrichment/
├── figures/
└── summary/
```

Interpretation:

* Final reliability requires both statistical and biological evidence.
* A reliable model should discriminate cases from controls, recover acceptable threshold-dependent performance, and produce interpretable probability estimates.
* Feature importance should identify signals that are not unique to a single learner unless a model-specific interpretation is intended.
* Reactome enrichment evaluates whether the final feature-importance structure maps onto biologically coherent pathways.

## Recommended execution order

Run the scripts in this order:

```r
source("scripts/01_hc_representation_diagnostics.R")
source("scripts/02_fsf_stability_conditioning.R")
source("scripts/03_benchmarking_ncv_ccv.R")
source("scripts/04_statistical_biological_reliability.R")
```

For large analyses, run each script in a clean R session.

## Reproducibility principles

TRACE-Omics follows these principles:

1. Use relative paths rather than local machine-specific paths.
2. Keep raw data outside version control.
3. Save intermediate objects required for downstream stages.
4. Store prediction-level outputs, not only scalar metrics.
5. Preserve `sessionInfo()` logs for every script.
6. Separate representation diagnostics from model evaluation.
7. Separate internal validation from cross-cohort validation.
8. Evaluate calibration and PR curves alongside AUROC.
9. Use final held-out test evaluation only after model selection.
10. Interpret final biomarkers using both model importance and pathway-level evidence.

## Public repository guidance

Do commit:

```text
README.md
LICENSE
CITATION.cff
.gitignore
config/
docs/
R/
scripts/
data/README.md
```

Do not commit:

```text
data/raw/
data/processed/
results/
*.rds
*.csv
*.png
*.tiff
*.pdf
*.html
```

## Minimal reproducible example

A minimal public example can be added later using synthetic data:

```text
data/example/
├── counts/
│   ├── cohort_A.csv
│   └── cohort_B.csv
└── metadata.rds
```

Synthetic data should preserve the expected structure but must not contain real patient-level values.

## Citation

If using TRACE-Omics, cite the repository and associated manuscript or thesis when available.

```
```

