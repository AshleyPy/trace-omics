# ======================================================================
# TB FULL-DATA DISEASE-SPECIFIC FINAL MODEL
# SEQUENTIAL | FAST | LEAKAGE-FREE | 80/20 FINAL ICE
# ======================================================================
#
# PURPOSE
#   1. Use the full TB population (expected n = 465).
#   2. Freeze one cohort x outcome-stratified 80/20 development/Test split.
#   3. Reproduce the equal-size representation diagnostics and Jaccard analysis.
#   4. Fast benchmark using ONLY VST + fSVA across:
#        FSFs: Information gain, RF importance, mRMR
#        Models: RBF SVM, Linear SVM, Random Forest
#        CBAs: None, cohort-stratified SMOTE, class weighting
#   5. Rank the benchmark primarily by BACC, then AUROC, robustness and Brier.
#   6. Carry all complete VST + fSVA benchmark candidates into strict CCV.
#   7. Lock ONE robust/generalizable pipeline before tuning.
#   8. Tune only that pipeline using nested CV, optimising BACC.
#   9. Quantify hyperparameter robustness and near-optimality.
#  10. Lock one robust hyperparameter configuration.
#  11. Fit on all development data and evaluate Final Ice exactly once.
#  12. Save a publication-ready trained model bundle and prediction helper.
#
# IMPORTANT
#   - Final Ice is never used for benchmark selection, CCV selection or tuning.
#   - Preprocessing, feature selection and SMOTE occur inside training data only.
#   - VST + fSVA is the only representation used throughout this definitive
#     disease-specific experiment, including diagnostics, benchmarking and CCV.
#   - Class weighting is learner-loss weighting, not resampling.
# ======================================================================

rm(list = ls(all.names = TRUE))
gc(verbose = FALSE)

required_packages <- c(
  "data.table", "DESeq2", "SummarizedExperiment", "sva", "matrixStats",
  "cluster", "neuroCombat", "edgeR", "mlr3verse", "mlr3filters",
  "mlr3learners", "mlr3extralearners", "ranger", "e1071", "pROC",
  "ggplot2", "RColorBrewer", "viridisLite", "scales", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop("Install these packages before running:\n", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
})

SEED <- 7832L
POSITIVE_CLASS <- "case"
NEGATIVE_CLASS <- "control"
TRAIN_FRACTION <- 0.80
TOP_N <- 50L
CV_FOLDS <- 5L
CV_REPEATS <- 2L
CCV_CANDIDATES <- 27L
MIN_TOTAL_TRAINING_COUNT <- 10L
MAX_SV <- 50L
MAX_PCA_K <- 50L
SILHOUETTE_NEAR_MAX <- 0.95
SMOTE_K <- 5L
TUNING_OUTER_FOLDS <- 5L
TUNING_INNER_FOLDS <- 3L
TUNING_EVALS <- 10L
FINAL_PERMUTATIONS <- 1000L
ROC_BOOTSTRAPS <- 500L

representation_ids <- c("VST_fSVA")
representation_labels <- c(
  VST_fSVA = "VST + fSVA"
)

fsf_ids <- c("information_gain", "cmim", "mrmr")
fsf_labels <- c(
  cmim = "CMIM",
  importance = "RF importance",
  mrmr = "mRMR"
)

learner_ids <- c("RBF_SVM", "Linear_SVM", "Random_Forest")
learner_labels <- c(
  RBF_SVM = "RBF SVM",
  Linear_SVM = "Linear SVM",
  Random_Forest = "Random Forest"
)

cba_ids <- c("None", "Cohort_SMOTE", "Class_Weighting")
cba_labels <- c(
  None = "None",
  Cohort_SMOTE = "Cohort-stratified SMOTE",
  Class_Weighting = "Class weighting"
)

set.seed(SEED)

seed_for <- function(...) {
  values <- as.numeric(c(...))
  weights <- seq_along(values) * 1009
  candidate <- as.numeric(SEED) + sum(values * weights)
  as.integer(candidate %% 2000000000)
}

user_home <- Sys.getenv("USERPROFILE", unset = path.expand("~"))
raw_tb_dir_candidates <- c(
  file.path(user_home, "Downloads", "RawIBDdata"),
  file.path(user_home, "OneDrive", "Downloads", "RawIBDdata")
)
raw_tb_dir_hits <- raw_tb_dir_candidates[dir.exists(raw_tb_dir_candidates)]
if (length(raw_tb_dir_hits) == 0L) {
  stop("RawIBDdata directory not found. Checked:\n", paste(raw_tb_dir_candidates, collapse = "\n"))
}
raw_tb_dir <- normalizePath(raw_tb_dir_hits[[1L]], winslash = "/", mustWork = TRUE)
counts_path <- file.path(raw_tb_dir, "cd_joined_numeric_matrix.rds")
metadata_path <- file.path(raw_tb_dir, "cd_metadata_495.rds")
if (!file.exists(counts_path)) stop("T1D raw-count file not found:\n", counts_path)
if (!file.exists(metadata_path)) stop("T1D metadata file not found:\n", metadata_path)

output_root <- file.path(user_home, "Desktop", "CD_FULL_DATA_FINAL_MODEL")
if (dir.exists(output_root)) {
  stop("Output directory already exists:\n", output_root,
       "\nRename or move the old directory before a clean rerun.")
}

audit_dir <- file.path(output_root, "00_audits")
split_dir <- file.path(output_root, "01_split")
diagnostic_dir <- file.path(output_root, "02_diagnostics")
concordance_dir <- file.path(output_root, "03_feature_concordance")
benchmark_dir <- file.path(output_root, "04_benchmark")
ccv_dir <- file.path(output_root, "05_CCV")
tuning_dir <- file.path(output_root, "06_tuning")
figure_dir <- file.path(output_root, "07_MAIN_FIGURES")
permutation_dir <- file.path(output_root, "08_PERMUTATION_TESTING")
final_ice_dir <- file.path(output_root, "99_FINAL_ICE")
bundle_dir <- file.path(output_root, "100_MODEL_BUNDLE")
for (d in c(output_root, audit_dir, split_dir, diagnostic_dir, concordance_dir,
            benchmark_dir, ccv_dir, tuning_dir, figure_dir, permutation_dir,
            final_ice_dir, bundle_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
