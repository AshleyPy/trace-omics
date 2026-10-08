# Dissertation methods-to-code map

This file maps the final dissertation Methods to the surviving TRACE-Omics code base. The mapping distinguishes the final statistical design from earlier development branches preserved for provenance.

| Dissertation analysis | Final design | Primary code location | Notes |
|---|---|---|---|
| Data harmonisation and partitioning | Disease-specific harmonisation; frozen 70:30 or 80:20 manifests depending on experiment | `scripts/01_representation_specification/` and `R/utils_io.R` | Exact historical manifests are external analytical state. |
| Preprocessing diagnostics | VST, VST + fSVA, VST + ComBat-Ref, ComBat-Seq + log1p with PCA/silhouette diagnostics | Stage 01 | Historical code uses `representation` as an internal variable name. |
| Feature-set characterisation | Jaccard agreement, VIF and effective rank; eight initial selectors, five retained | Stage 02 | Matches the feature-set analyses in Appendix C. |
| Principal 105-pipeline benchmark | 3 preprocessing strategies × 5 FS methods × 7 classifiers; 5 outer folds; 10 inner folds when tuned | Stage 03 functions plus retained historical benchmark code under `legacy/` | Exact original fold/tuning archives are required for numerical identity. |
| Equal-size sensitivity | TB/RA/UC/CD; 100 repeated study-aware draws; n=82 with 41 cases/41 controls | Final specification in `config/final_dissertation.yml`; related historical code under `legacy/equal_size/` | Older 66/16 ICE branches are provenance only and are not the final sensitivity estimand. |
| Study transportability | LOSO-CV within 80% development; VST or log2(CPM+1); RF importance; 25/50/100 genes; weighted/unweighted RF | `scripts/05_study_transportability.R`, `R/study_models.R`, Stage 03 Experiment 3 | LOSO-CV configuration selection and target-specific evaluation share the available study-exclusion design; this limitation is reported in Chapter 5. |
| Matched LightGBM–TabPFN | Same observations and 100 RF-selected genes; 5 outer/3 inner folds; threshold 0.50 | `scripts/08_matched_architectures.R`, Stage 03/04 Experiment 4, `legacy/matched_architectures/` | Requires the recorded Python environment and TabPFN checkpoint. |
| Pooled Models A/B | Fixed 4,450 genes; Model A n=1,413; Model B n=1,736; same 762-observation holdout | Stage 03 Experiment 2 and Stage 04 Experiment 2 | Sample size and disease composition change together. |
| Pooled analytical sensitivity | VST + ComBat-Ref, Boruta 50, observed distribution, LightGBM | `legacy/pooled_combatref/` plus Stage 04 utilities | Preprocessing and FS change together. |
| Complete disease exclusion | LODO-CV; fixed primary pooled system: VST + fSVA, mRMR 50, observed distribution, LightGBM | Stage 04 Experiment 2 functions | fSVA is carried forward because the analytical system is held fixed while the disease boundary changes. |
| Predictive uncertainty | DeLong AUROC CI; paired DeLong for matched predictions; stratified bootstrap for reported non-AUROC metrics | `R/core.R`, Stage 04 | Fold/study/disease SD describes resampling-unit variability rather than an integrated CI. |
| Molecular interpretation | Permutation importance, rank/Jaccard agreement, ALE, model restoration audit | Stage 04 and `legacy/matched_architectures/` | Predictive reliance is kept separate from biological validation. |

## Terminology

The final dissertation uses **preprocessing strategy**, **feature selection (FS)**, **class imbalance correction (CIC)**, **population generalisability**, **study transportability**, **disease transportability** and **sensitivity analysis**. Internal code variable names inherited from development, especially `representation`, are preserved where renaming would change provenance-bearing source.
