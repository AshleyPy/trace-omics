# RELEASE NOTE -----------------------------------------------------------------
# This source is a recovered comprehensive TRACE-Omics implementation from the
# dissertation code-development history. The final statistical design is defined
# in config/final_dissertation.yml and docs/METHODS_TO_CODE.md. Earlier branches
# retained inside this file are provenance-bearing and may differ from the final
# reported sensitivity estimand. See docs/REPRODUCTION_STATUS.md before rerunning.
# -----------------------------------------------------------------------------
# ==============================================================================
# TRACE-Omics
# Stage 03: Model development and optimisation
# ==============================================================================
#
# File:
#   scripts/03_model_development_optimisation/03_model_development_optimisation.R
#
# Framework:
#   Transcriptomic Robustness Assessment and Cross-cohort Evaluation for
#   Omics learning (TRACE-Omics)
#
# Purpose:
#   Execute the complete development-only modelling and optimisation stage used
#   across the four controlled TRACE-Omics experiments.
#
#   This stage consumes the canonical data/manifests created by Stage 01 and the
#   feature-selection machinery defined in Stage 02. It produces development
#   predictions, optimisation archives and locked analytical specifications for
#   Stage 04.
#
#   Experiment 1 - sample-support perturbation
#     * 10-fold x 2-repeat cohort/outcome-stratified CV.
#     * Fold-internal representation fitting and feature selection.
#     * Five feature selectors x eight fixed learners across three representations.
#     * Development-only representation-sensitivity reference:
#         RF importance | 50 genes | LightGBM.
#
#   Experiment 2 - disease-composition perturbation
#     * Fixed 4,450-gene cross-disease space.
#     * Model A and Model B fitted independently.
#     * 2 feature selectors x 2 class-handling states x 4 learner families.
#     * 5 outer folds, 10 inner folds, 20 MBO evaluations.
#     * Balanced accuracy used as the tuning objective.
#     * Preprocessing, feature selection and SMOTE repeated inside every fitting
#       partition.
#     * Locked analytical pipeline:
#         VST + fSVA | mRMR | None | LightGBM.
#     * A final development-only hyperparameter search and OOF-derived Youden
#       threshold are saved for subsequent Stage-04 fitting.
#
#   Experiment 3 - cohort-composition perturbation
#     * Development-only LOCO optimisation.
#     * VST only versus log2(CPM + 1).
#     * Training-restricted RF permutation feature ranking.
#     * RBF-SVM four-configuration screen and Random Forest screen.
#     * 25 / 50 / 100 genes x unweighted / class-weighted models.
#     * Final Random Forest configuration ranked by:
#         1. highest worst-cohort AUROC;
#         2. highest mean LOCO AUROC;
#         3. lowest between-cohort AUROC SD;
#         4. highest mean LOCO BACC.
#     * Development LOCO threshold selected from 0.05-0.95 by:
#         1. highest mean cohort BACC;
#         2. highest worst-cohort BACC;
#         3. closest threshold to 0.50.
#
#   Experiment 4 - model-architecture perturbation
#     * Matched log2(CPM + 1) + RF-importance top-50 inputs.
#     * 5 outer folds, 3 inner folds.
#     * Prespecified ten-configuration LightGBM grid.
#     * Default pretrained TabPFN v0.1.11.
#     * LightGBM and TabPFN receive the same outer-training observations,
#       same outer-validation observations, same 50 genes and same gene order.
#     * Final LightGBM configuration is locked from development resampling using:
#         mean inner AUROC -> lower SD -> stronger worst-fold inner AUROC.
#
# Stage boundary:
#   Stage 03 never evaluates the frozen final held-out population. It prepares
#   and locks all model-development choices. Stage 04 performs final frozen
#   held-out evaluation, generalisation summaries, calibration, statistical
#   comparisons and molecular reliability analyses.
#
# Example:
#
#   Rscript scripts/03_model_development_optimisation/03_model_development_optimisation.R
#
#   Rscript scripts/03_model_development_optimisation/03_model_development_optimisation.R \
#     --config=config/example_config.yml \
#     --mode=all \
#     --diseases=TB,RA,UC,CD,T1D \
#     --strict=false
#
# ==============================================================================


# ==============================================================================
# 0. COMMAND-LINE OPTIONS
# ==============================================================================

parse_named_arguments <- function(x) {
  output <- list()

  for (arg in x) {
    if (!grepl("^--", arg)) {
      next
    }

    pieces <- strsplit(
      sub("^--", "", arg),
      "=",
      fixed = TRUE
    )[[1]]

    key <- pieces[[1]]
    value <- if (length(pieces) >= 2L) {
      paste(pieces[-1], collapse = "=")
    } else {
      "true"
    }

    output[[key]] <- value
  }

  output
}

parse_boolean <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }

  tolower(as.character(x)) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}

args <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE)
)

CONFIG_PATH <- if (!is.null(args$config)) {
  args$config
} else {
  "config/example_config.yml"
}

MODE_TEXT <- if (!is.null(args$mode)) {
  args$mode
} else {
  "all"
}

DISEASE_TEXT <- if (!is.null(args$diseases)) {
  args$diseases
} else {
  "TB,RA,UC,CD,T1D"
}

STRICT_AUDITS <- parse_boolean(
  args$strict,
  default = FALSE
)

RUN_MODES <- unique(
  trimws(
    strsplit(
      MODE_TEXT,
      ",",
      fixed = TRUE
    )[[1]]
  )
)

REQUESTED_DISEASES <- unique(
  trimws(
    strsplit(
      DISEASE_TEXT,
      ",",
      fixed = TRUE
    )[[1]]
  )
)

VALID_MODES <- c(
  "all",
  "experiment1",
  "experiment2",
  "experiment3",
  "experiment4"
)

if (!all(RUN_MODES %in% VALID_MODES)) {
  stop(
    "Unknown --mode value(s): ",
    paste(
      setdiff(
        RUN_MODES,
        VALID_MODES
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

run_requested <- function(mode) {
  "all" %in% RUN_MODES ||
    mode %in% RUN_MODES
}


# ==============================================================================
# 1. CORE UTILITIES, PACKAGES AND CONFIGURATION
# ==============================================================================

UTILS_PATH <- file.path(
  "R",
  "utils_io.R"
)

if (!file.exists(UTILS_PATH)) {
  stop(
    "TRACE-Omics utility file not found: ",
    UTILS_PATH,
    "\nRun this script from the repository root.",
    call. = FALSE
  )
}

source(
  UTILS_PATH
)

require_packages(
  c(
    "yaml",
    "data.table",
    "matrixStats",
    "mlr3",
    "mlr3learners",
    "mlr3extralearners",
    "paradox",
    "bbotk",
    "mlr3mbo",
    "ranger",
    "e1071",
    "lightgbm",
    "xgboost",
    "C50",
    "naivebayes",
    "FNN",
    "pROC",
    "PRROC"
  ),
  bioc = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
  library(mlr3)
  library(paradox)
  library(bbotk)
  library(mlr3mbo)
})

config <- read_trace_config(
  CONFIG_PATH
)

GLOBAL_SEED <- as.integer(
  get_config_value(
    config,
    "project.seed",
    7832
  )
)

POSITIVE_CLASS <- as.character(
  get_config_value(
    config,
    "project.positive_class",
    "case"
  )
)

NEGATIVE_CLASS <- as.character(
  get_config_value(
    config,
    "project.negative_class",
    "control"
  )
)

ALL_DISEASES <- c(
  "TB",
  "RA",
  "UC",
  "CD",
  "T1D"
)

if (!all(REQUESTED_DISEASES %in% ALL_DISEASES)) {
  stop(
    "Unknown disease identifier(s): ",
    paste(
      setdiff(
        REQUESTED_DISEASES,
        ALL_DISEASES
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

set_trace_seed(
  GLOBAL_SEED
)


# ==============================================================================
# 2. IMPORT FUNCTION DEFINITIONS FROM STAGES 01 AND 02
#
# Only the function-definition prefixes are imported. The execution blocks at
# the bottom of the prior-stage scripts are deliberately excluded.
# ==============================================================================

STAGE1_SCRIPT <- file.path(
  "scripts",
  "01_representation_specification",
  "01_representation_specification.R"
)

STAGE2_SCRIPT <- file.path(
  "scripts",
  "02_feature_set_characterisation",
  "02_feature_set_characterisation.R"
)

source_prefix_before_marker <- function(path,
                                        marker,
                                        environment) {
  check_file_exists(
    path,
    paste0(
      "TRACE-Omics prerequisite script: ",
      basename(path)
    )
  )

  lines <- readLines(
    path,
    warn = FALSE,
    encoding = "UTF-8"
  )

  marker_line <- grep(
    marker,
    lines,
    fixed = TRUE
  )

  if (length(marker_line) == 0L) {
    stop(
      "Marker not found in ",
      path,
      ": ",
      marker,
      call. = FALSE
    )
  }

  keep_to <- marker_line[[1]] - 1L

  temporary_file <- tempfile(
    fileext = ".R"
  )

  on.exit(
    unlink(temporary_file),
    add = TRUE
  )

  writeLines(
    lines[seq_len(keep_to)],
    temporary_file,
    useBytes = TRUE
  )

  sys.source(
    temporary_file,
    envir = environment
  )

  invisible(environment)
}

stage1 <- new.env(
  parent = globalenv()
)

stage2 <- new.env(
  parent = globalenv()
)

source_prefix_before_marker(
  STAGE1_SCRIPT,
  "# 20. LOAD ALL REQUESTED DISEASES AND CREATE FROZEN GENERAL MANIFESTS",
  stage1
)

source_prefix_before_marker(
  STAGE2_SCRIPT,
  "# 11. FULL-DEVELOPMENT FEATURE-SET CHARACTERISATION",
  stage2
)

required_stage1_functions <- c(
  "fit_apply_representation",
  "fit_apply_vst_fsva",
  "fit_apply_frozen_vst",
  "fit_apply_log2_cpm_plus_1"
)

required_stage2_functions <- c(
  "select_ranked_features",
  "fit_apply_log2_cpm_feature_stage"
)

for (function_name in required_stage1_functions) {
  if (!exists(
    function_name,
    envir = stage1,
    inherits = FALSE
  )) {
    stop(
      "Stage 01 function unavailable after import: ",
      function_name,
      call. = FALSE
    )
  }
}

for (function_name in required_stage2_functions) {
  if (!exists(
    function_name,
    envir = stage2,
    inherits = FALSE
  )) {
    stop(
      "Stage 02 function unavailable after import: ",
      function_name,
      call. = FALSE
    )
  }
}


# ==============================================================================
# 3. OUTPUT ROOTS AND REPRODUCIBILITY RECORDS
# ==============================================================================

stage3_example <- get_config_value(
  config,
  "paths.stage_results.model_development_optimisation",
  "results/03_model_development_optimisation/TB"
)

STAGE3_ROOT <- dirname(
  stage3_example
)

AUDIT_ROOT <- file.path(
  STAGE3_ROOT,
  "00_audits"
)

SESSION_ROOT <- file.path(
  STAGE3_ROOT,
  "99_session_info"
)

make_dir(
  STAGE3_ROOT
)

make_dir(
  AUDIT_ROOT
)

make_dir(
  SESSION_ROOT
)

stage_source_manifest <- data.table(
  stage = c(
    "Stage01",
    "Stage02"
  ),
  path = c(
    STAGE1_SCRIPT,
    STAGE2_SCRIPT
  ),
  md5 = unname(
    tools::md5sum(
      c(
        STAGE1_SCRIPT,
        STAGE2_SCRIPT
      )
    )
  )
)

safe_write_csv(
  stage_source_manifest,
  file.path(
    AUDIT_ROOT,
    "TRACE_Omics_stage03_prerequisite_md5.csv"
  )
)

message_section(
  "TRACE-Omics Stage 03: Model development and optimisation"
)

message(
  "Configuration: ",
  CONFIG_PATH
)

message(
  "Modes:         ",
  paste(
    RUN_MODES,
    collapse = ", "
  )
)

message(
  "Diseases:      ",
  paste(
    REQUESTED_DISEASES,
    collapse = ", "
  )
)

message(
  "Global seed:   ",
  GLOBAL_SEED
)


# ==============================================================================
# 4. CANONICAL DATA AND MANIFEST LOADERS
# ==============================================================================

load_canonical_disease <- function(disease_id) {
  base_dir <- file.path(
    "data",
    "processed",
    disease_id,
    "representation_specification"
  )

  count_path <- file.path(
    base_dir,
    paste0(
      disease_id,
      "_canonical_raw_counts.rds"
    )
  )

  metadata_path <- file.path(
    base_dir,
    paste0(
      disease_id,
      "_canonical_metadata.rds"
    )
  )

  manifest_path <- file.path(
    "data",
    "processed",
    disease_id,
    "manifests",
    paste0(
      disease_id,
      "_development_heldout_manifest.csv"
    )
  )

  check_file_exists(
    count_path,
    paste0(
      disease_id,
      " canonical raw counts"
    )
  )

  check_file_exists(
    metadata_path,
    paste0(
      disease_id,
      " canonical metadata"
    )
  )

  check_file_exists(
    manifest_path,
    paste0(
      disease_id,
      " development/held-out manifest"
    )
  )

  counts <- coerce_numeric_matrix(
    safe_read_rds(
      count_path
    ),
    paste0(
      disease_id,
      " canonical counts"
    )
  )

  metadata <- as.data.frame(
    safe_read_rds(
      metadata_path
    ),
    check.names = FALSE
  )

  manifest <- as.data.frame(
    safe_read_csv(
      manifest_path
    ),
    check.names = FALSE
  )

  check_required_columns(
    metadata,
    c(
      "sample_id",
      "cohort_group",
      "outcome"
    ),
    paste0(
      disease_id,
      " canonical metadata"
    )
  )

  check_required_columns(
    manifest,
    c(
      "sample_id",
      "role"
    ),
    paste0(
      disease_id,
      " partition manifest"
    )
  )

  metadata$sample_id <- as.character(
    metadata$sample_id
  )

  metadata$cohort_group <- factor(
    as.character(
      metadata$cohort_group
    )
  )

  metadata$outcome <- factor(
    as.character(
      metadata$outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  development_ids <- manifest$sample_id[
    manifest$role == "development"
  ]

  development_metadata <- metadata[
    match(
      development_ids,
      metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(
    development_metadata$sample_id
  )) {
    stop(
      disease_id,
      " development manifest contains sample IDs absent from canonical metadata.",
      call. = FALSE
    )
  }

  development_counts <- counts[
    ,
    development_metadata$sample_id,
    drop = FALSE
  ]

  list(
    disease = disease_id,
    counts = counts,
    metadata = metadata,
    manifest = manifest,
    development_counts = development_counts,
    development_metadata = development_metadata
  )
}

load_equal_size_development <- function(disease_id) {
  manifest_path <- file.path(
    dirname(
      get_config_value(
        config,
        "paths.stage_results.representation_specification",
        "results/01_representation_specification/TB"
      )
    ),
    "02_experiment_1_sample_support",
    disease_id,
    "00_manifests",
    paste0(
      disease_id,
      "_equal_size_82_manifest.csv"
    )
  )

  check_file_exists(
    manifest_path,
    paste0(
      disease_id,
      " equal-size manifest"
    )
  )

  canonical <- load_canonical_disease(
    disease_id
  )

  manifest <- as.data.frame(
    safe_read_csv(
      manifest_path
    ),
    check.names = FALSE
  )

  development_ids <- manifest$sample_id[
    manifest$role == "development"
  ]

  metadata <- canonical$metadata[
    match(
      development_ids,
      canonical$metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  counts <- canonical$counts[
    ,
    metadata$sample_id,
    drop = FALSE
  ]

  if (
    nrow(metadata) != 66L ||
      sum(
        metadata$outcome ==
          POSITIVE_CLASS
      ) != 33L ||
      sum(
        metadata$outcome ==
          NEGATIVE_CLASS
      ) != 33L
  ) {
    stop(
      disease_id,
      " equal-size development population differs from the prespecified 66 (33/33) design.",
      call. = FALSE
    )
  }

  list(
    counts = counts,
    metadata = metadata
  )
}


# ==============================================================================
# 5. FOLD CONSTRUCTION
# ==============================================================================

make_joint_stratified_folds <- function(metadata,
                                        number_of_folds,
                                        seed,
                                        include_disease = FALSE) {
  if (number_of_folds < 2L) {
    stop(
      "number_of_folds must be >= 2.",
      call. = FALSE
    )
  }

  if (isTRUE(include_disease)) {
    strata <- interaction(
      metadata$disease,
      metadata$cohort_group,
      metadata$outcome,
      drop = TRUE,
      lex.order = TRUE
    )
  } else {
    strata <- interaction(
      metadata$cohort_group,
      metadata$outcome,
      drop = TRUE,
      lex.order = TRUE
    )
  }

  fold_id <- integer(
    nrow(metadata)
  )

  stratum_levels <- unique(
    as.character(
      strata
    )
  )

  for (stratum_index in seq_along(stratum_levels)) {
    stratum_name <- stratum_levels[[stratum_index]]

    rows <- which(
      as.character(strata) ==
        stratum_name
    )

    set_trace_subseed(
      seed,
      "fold",
      number_of_folds,
      stratum_name
    )

    rows <- sample(
      rows,
      length(rows),
      replace = FALSE
    )

    assigned <- rep(
      seq_len(
        number_of_folds
      ),
      length.out = length(rows)
    )

    fold_id[rows] <- assigned
  }

  fold_id
}

make_repeated_fold_map <- function(metadata,
                                   folds,
                                   repeats,
                                   seed,
                                   include_disease = FALSE) {
  rows <- list()

  for (repeat_index in seq_len(repeats)) {
    current_fold <- make_joint_stratified_folds(
      metadata = metadata,
      number_of_folds = folds,
      seed = make_trace_subseed(
        seed,
        "repeat",
        `repeat` = repeat_index
      ),
      include_disease = include_disease
    )

    rows[[repeat_index]] <- data.table(
      sample_id = metadata$`sample_id`,
      `repeat` = repeat_index,
      fold = current_fold,
      cohort = as.character(
        metadata$cohort_group
      ),
      outcome = as.character(
        metadata$outcome
      )
    )

    if ("disease" %in% names(metadata)) {
      rows[[repeat_index]][
        ,
        disease := as.character(
          metadata$disease
        )
      ]
    }
  }

  rbindlist(
    rows,
    fill = TRUE
  )
}


# ==============================================================================
# 6. PERFORMANCE METRICS
# ==============================================================================

safe_auroc <- function(truth,
                       probability) {
  truth <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  valid <- (
    !is.na(truth) &
      is.finite(probability)
  )

  truth <- truth[valid]
  probability <- probability[valid]

  if (
    length(unique(truth)) < 2L ||
      length(probability) < 2L
  ) {
    return(NA_real_)
  }

  roc_object <- tryCatch(
    pROC::roc(
      response = truth,
      predictor = probability,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      ),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(roc_object)) {
    return(NA_real_)
  }

  as.numeric(
    pROC::auc(
      roc_object
    )
  )
}

safe_auprc <- function(truth,
                       probability) {
  truth <- as.character(
    truth
  )

  valid <- (
    truth %in% c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ) &
      is.finite(probability)
  )

  truth <- truth[valid]
  probability <- probability[valid]

  positive_scores <- probability[
    truth ==
      POSITIVE_CLASS
  ]

  negative_scores <- probability[
    truth ==
      NEGATIVE_CLASS
  ]

  if (
    length(positive_scores) == 0L ||
      length(negative_scores) == 0L
  ) {
    return(NA_real_)
  }

  curve <- tryCatch(
    PRROC::pr.curve(
      scores.class0 = positive_scores,
      scores.class1 = negative_scores,
      curve = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(curve)) {
    return(NA_real_)
  }

  as.numeric(
    curve$auc.integral
  )
}

classification_metrics <- function(truth,
                                   probability,
                                   threshold = 0.50) {
  truth <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  valid <- (
    !is.na(truth) &
      is.finite(probability)
  )

  truth <- truth[valid]
  probability <- probability[valid]

  if (length(probability) == 0L) {
    return(
      data.table(
        n = 0L,
        prevalence = NA_real_,
        AUROC = NA_real_,
        AUPRC = NA_real_,
        BACC = NA_real_,
        Sensitivity = NA_real_,
        Specificity = NA_real_,
        Brier = NA_real_
      )
    )
  }

  predicted <- factor(
    ifelse(
      probability >= threshold,
      POSITIVE_CLASS,
      NEGATIVE_CLASS
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  tp <- sum(
    predicted ==
      POSITIVE_CLASS &
      truth ==
      POSITIVE_CLASS
  )

  tn <- sum(
    predicted ==
      NEGATIVE_CLASS &
      truth ==
      NEGATIVE_CLASS
  )

  fp <- sum(
    predicted ==
      POSITIVE_CLASS &
      truth ==
      NEGATIVE_CLASS
  )

  fn <- sum(
    predicted ==
      NEGATIVE_CLASS &
      truth ==
      POSITIVE_CLASS
  )

  sensitivity <- if (
    tp + fn > 0L
  ) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  specificity <- if (
    tn + fp > 0L
  ) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  bacc <- mean(
    c(
      sensitivity,
      specificity
    ),
    na.rm = TRUE
  )

  outcome_numeric <- as.integer(
    truth ==
      POSITIVE_CLASS
  )

  data.table(
    n = length(truth),
    prevalence = mean(
      outcome_numeric
    ),
    AUROC = safe_auroc(
      truth,
      probability
    ),
    AUPRC = safe_auprc(
      truth,
      probability
    ),
    BACC = bacc,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Brier = mean(
      (
        probability -
          outcome_numeric
      )^2
    )
  )
}


# ==============================================================================
# 7. CLASSIFICATION TASK AND FIXED LEARNER HELPERS
# ==============================================================================

make_mlr3_task <- function(expression_matrix,
                           metadata,
                           task_id) {
  expression_matrix <- coerce_numeric_matrix(
    expression_matrix,
    paste0(
      task_id,
      " expression matrix"
    )
  )

  sample_ids <- intersect(
    colnames(expression_matrix),
    metadata$sample_id
  )

  expression_matrix <- expression_matrix[
    ,
    sample_ids,
    drop = FALSE
  ]

  metadata <- metadata[
    match(
      sample_ids,
      metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  original_genes <- rownames(
    expression_matrix
  )

  safe_genes <- make.unique(
    make.names(
      original_genes,
      unique = TRUE
    )
  )

  x <- as.data.table(
    t(
      expression_matrix
    )
  )

  setnames(
    x,
    safe_genes
  )

  x[
    ,
    outcome := factor(
      as.character(
        metadata$outcome
      ),
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    )
  ]

  task <- mlr3::TaskClassif$new(
    id = task_id,
    backend = x,
    target = "outcome",
    positive = POSITIVE_CLASS
  )

  list(
    task = task,
    map = data.table(
      original_gene = original_genes,
      safe_gene = safe_genes
    ),
    metadata = metadata
  )
}

make_fixed_learner <- function(learner_id,
                               seed) {
  if (learner_id == "Naive_Bayes") {
    return(
      mlr3::lrn(
        "classif.naive_bayes",
        predict_type = "prob",
        laplace = 0
      )
    )
  }

  if (learner_id == "Linear_SVM") {
    return(
      mlr3::lrn(
        "classif.svm",
        predict_type = "prob",
        type = "C-classification",
        kernel = "linear",
        cost = 1,
        scale = TRUE,
        cachesize = 500,
        tolerance = 1e-3,
        shrinking = TRUE
      )
    )
  }

  if (learner_id == "RBF_SVM") {
    return(
      mlr3::lrn(
        "classif.svm",
        predict_type = "prob",
        type = "C-classification",
        kernel = "radial",
        scale = TRUE,
        cachesize = 500,
        tolerance = 1e-3,
        shrinking = TRUE
      )
    )
  }

  if (learner_id == "C5_0") {
    return(
      mlr3::lrn(
        "classif.C50",
        predict_type = "prob",
        trials = 1L,
        rules = FALSE,
        winnow = FALSE
      )
    )
  }

  if (learner_id == "Random_Forest") {
    return(
      mlr3::lrn(
        "classif.ranger",
        predict_type = "prob",
        num.trees = 500L,
        num.threads = 1L,
        importance = "none",
        save.memory = TRUE,
        seed = as.integer(seed)
      )
    )
  }

  if (learner_id == "LightGBM") {
    return(
      mlr3::lrn(
        "classif.lightgbm",
        predict_type = "prob",
        num_iterations = 250L,
        min_data_in_leaf = 5L,
        bagging_fraction = 0.80,
        bagging_freq = 1L,
        num_threads = 1L,
        device_type = "cpu",
        deterministic = TRUE,
        seed = as.integer(seed)
      )
    )
  }

  if (learner_id == "XGBoost") {
    return(
      mlr3::lrn(
        "classif.xgboost",
        predict_type = "prob",
        nrounds = 250L,
        min_child_weight = 1,
        colsample_bytree = 0.80,
        subsample = 0.80,
        nthread = 1L,
        eval_metric = "logloss",
        seed = as.integer(seed)
      )
    )
  }

  stop(
    "Unknown fixed learner: ",
    learner_id,
    call. = FALSE
  )
}

fit_predict_fixed_learner <- function(learner_id,
                                      train_matrix,
                                      train_metadata,
                                      validation_matrix,
                                      validation_metadata,
                                      seed) {
  if (learner_id == "Featureless") {
    prevalence <- mean(
      train_metadata$outcome ==
        POSITIVE_CLASS
    )

    return(
      rep(
        prevalence,
        nrow(validation_metadata)
      )
    )
  }

  training_task_info <- make_mlr3_task(
    expression_matrix = train_matrix,
    metadata = train_metadata,
    task_id = paste0(
      "train_",
      safe_file_name(
        learner_id
      )
    )
  )

  validation_task_info <- make_mlr3_task(
    expression_matrix = validation_matrix,
    metadata = validation_metadata,
    task_id = paste0(
      "valid_",
      safe_file_name(
        learner_id
      )
    )
  )

  learner <- make_fixed_learner(
    learner_id = learner_id,
    seed = seed
  )

  set.seed(
    seed
  )

  learner$train(
    training_task_info$task
  )

  prediction <- learner$predict(
    validation_task_info$task
  )

  probability <- prediction$prob[
    ,
    POSITIVE_CLASS
  ]

  as.numeric(
    probability
  )
}


# ==============================================================================
# 8. COHORT/DISEASE-STRATIFIED SMOTE
# ==============================================================================

smote_one_stratum <- function(expression_matrix,
                              metadata,
                              k = 5L,
                              seed = GLOBAL_SEED) {
  class_counts <- table(
    metadata$outcome
  )

  present <- names(
    class_counts[
      class_counts > 0L
    ]
  )

  if (length(present) < 2L) {
    return(
      list(
        matrix = expression_matrix,
        metadata = metadata,
        synthetic_n = 0L,
        note = "Single outcome class in stratum"
      )
    )
  }

  minority_class <- names(
    which.min(
      class_counts[
        class_counts > 0L
      ]
    )
  )

  majority_class <- names(
    which.max(
      class_counts[
        class_counts > 0L
      ]
    )
  )

  minority_ids <- metadata$sample_id[
    metadata$outcome ==
      minority_class
  ]

  minority_n <- length(
    minority_ids
  )

  majority_n <- sum(
    metadata$outcome ==
      majority_class
  )

  synthetic_n <- majority_n -
    minority_n

  if (synthetic_n <= 0L) {
    return(
      list(
        matrix = expression_matrix,
        metadata = metadata,
        synthetic_n = 0L,
        note = NA_character_
      )
    )
  }

  x_minority <- t(
    expression_matrix[
      ,
      minority_ids,
      drop = FALSE
    ]
  )

  storage.mode(
    x_minority
  ) <- "double"

  set.seed(
    seed
  )

  if (minority_n < 2L) {
    synthetic_matrix <- matrix(
      rep(
        x_minority[
          1L,
          ,
          drop = TRUE
        ],
        synthetic_n
      ),
      nrow = synthetic_n,
      byrow = TRUE
    )

    colnames(
      synthetic_matrix
    ) <- colnames(
      x_minority
    )
  } else {
    effective_k <- min(
      as.integer(k),
      minority_n - 1L
    )

    neighbour_indices <- FNN::get.knn(
      x_minority,
      k = effective_k
    )$nn.index

    synthetic_matrix <- matrix(
      NA_real_,
      nrow = synthetic_n,
      ncol = ncol(x_minority)
    )

    colnames(
      synthetic_matrix
    ) <- colnames(
      x_minority
    )

    for (synthetic_index in seq_len(synthetic_n)) {
      seed_row <- sample.int(
        minority_n,
        size = 1L
      )

      neighbour_row <- sample(
        neighbour_indices[
          seed_row,

        ],
        size = 1L
      )

      interpolation <- stats::runif(
        1L
      )

      synthetic_matrix[
        synthetic_index,

      ] <- x_minority[
        seed_row,

      ] +
        interpolation *
          (
            x_minority[
              neighbour_row,

            ] -
              x_minority[
                seed_row,

              ]
          )
    }
  }

  synthetic_ids <- paste0(
    "SMOTE_",
    make_trace_subseed(
      seed,
      "synthetic"
    ),
    "_",
    seq_len(synthetic_n)
  )

  rownames(
    synthetic_matrix
  ) <- synthetic_ids

  synthetic_metadata <- metadata[
    rep(
      which(
        metadata$outcome ==
          minority_class
      )[[1]],
      synthetic_n
    ),
    ,
    drop = FALSE
  ]

  synthetic_metadata$sample_id <- synthetic_ids
  synthetic_metadata$outcome <- factor(
    rep(
      minority_class,
      synthetic_n
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  combined_matrix <- cbind(
    expression_matrix,
    t(
      synthetic_matrix
    )
  )

  combined_metadata <- rbind(
    metadata,
    synthetic_metadata
  )

  combined_metadata <- combined_metadata[
    match(
      colnames(combined_matrix),
      combined_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  list(
    matrix = combined_matrix,
    metadata = combined_metadata,
    synthetic_n = synthetic_n,
    note = if (minority_n < 2L) {
      "Minority duplication used because one minority observation was available"
    } else {
      NA_character_
    }
  )
}

stratified_smote <- function(expression_matrix,
                             metadata,
                             k = 5L,
                             seed = GLOBAL_SEED,
                             include_disease = TRUE) {
  metadata <- as.data.frame(
    metadata,
    check.names = FALSE
  )

  if (
    isTRUE(include_disease) &&
      "disease" %in% colnames(metadata)
  ) {
    strata <- interaction(
      metadata$disease,
      metadata$cohort_group,
      drop = TRUE,
      lex.order = TRUE
    )
  } else {
    strata <- factor(
      as.character(
        metadata$cohort_group
      )
    )
  }

  matrices <- list()
  metadata_parts <- list()
  audit_rows <- list()

  stratum_levels <- unique(
    as.character(
      strata
    )
  )

  for (stratum_index in seq_along(stratum_levels)) {
    stratum_name <- stratum_levels[[stratum_index]]

    rows <- which(
      as.character(strata) ==
        stratum_name
    )

    current_metadata <- metadata[
      rows,
      ,
      drop = FALSE
    ]

    current_matrix <- expression_matrix[
      ,
      current_metadata$sample_id,
      drop = FALSE
    ]

    result <- smote_one_stratum(
      expression_matrix = current_matrix,
      metadata = current_metadata,
      k = k,
      seed = make_trace_subseed(
        seed,
        stratum_name
      )
    )

    matrices[[stratum_name]] <- result$matrix
    metadata_parts[[stratum_name]] <- result$metadata

    audit_rows[[stratum_name]] <- data.table(
      stratum = stratum_name,
      original_n = nrow(
        current_metadata
      ),
      synthetic_n = result$synthetic_n,
      final_n = nrow(
        result$metadata
      ),
      note = result$note
    )
  }

  combined_matrix <- do.call(
    cbind,
    matrices
  )

  combined_metadata <- do.call(
    rbind,
    metadata_parts
  )

  combined_metadata <- as.data.frame(
    combined_metadata,
    check.names = FALSE
  )

  combined_metadata <- combined_metadata[
    match(
      colnames(combined_matrix),
      combined_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  set_trace_subseed(
    seed,
    "smote_shuffle"
  )

  shuffled <- sample(
    seq_len(
      nrow(combined_metadata)
    ),
    replace = FALSE
  )

  combined_metadata <- combined_metadata[
    shuffled,
    ,
    drop = FALSE
  ]

  combined_matrix <- combined_matrix[
    ,
    combined_metadata$sample_id,
    drop = FALSE
  ]

  list(
    matrix = combined_matrix,
    metadata = combined_metadata,
    audit = rbindlist(
      audit_rows,
      fill = TRUE
    )
  )
}


# ==============================================================================
# 9. TUNABLE LEARNERS AND MBO SEARCH SPACES
# ==============================================================================

CROSS_DISEASE_SEARCH_SPACES <- list(
  Random_Forest = paradox::ps(
    num.trees = paradox::p_int(
      lower = 200L,
      upper = 500L
    ),
    mtry.ratio = paradox::p_dbl(
      lower = 0.05,
      upper = 0.75
    ),
    min.node.size = paradox::p_int(
      lower = 1L,
      upper = 15L
    )
  ),
  LightGBM = paradox::ps(
    learning_rate = paradox::p_dbl(
      lower = 0.01,
      upper = 0.25
    ),
    num_leaves = paradox::p_int(
      lower = 8L,
      upper = 100L
    ),
    feature_fraction = paradox::p_dbl(
      lower = 0.50,
      upper = 1.00
    ),
    bagging_fraction = paradox::p_dbl(
      lower = 0.50,
      upper = 1.00
    ),
    min_data_in_leaf = paradox::p_int(
      lower = 5L,
      upper = 50L
    ),
    num_iterations = paradox::p_int(
      lower = 100L,
      upper = 500L
    )
  ),
  XGBoost = paradox::ps(
    eta = paradox::p_dbl(
      lower = 0.01,
      upper = 0.30
    ),
    max_depth = paradox::p_int(
      lower = 2L,
      upper = 8L
    ),
    subsample = paradox::p_dbl(
      lower = 0.60,
      upper = 1.00
    ),
    colsample_bytree = paradox::p_dbl(
      lower = 0.50,
      upper = 1.00
    ),
    nrounds = paradox::p_int(
      lower = 100L,
      upper = 500L
    )
  ),
  RBF_SVM = paradox::ps(
    cost = paradox::p_dbl(
      lower = 1e-2,
      upper = 1e2,
      logscale = TRUE
    ),
    gamma = paradox::p_dbl(
      lower = 1e-5,
      upper = 2,
      logscale = TRUE
    )
  )
)

make_tunable_learner <- function(learner_id,
                                 parameters,
                                 seed,
                                 class_weights = NULL) {
  if (learner_id == "Random_Forest") {
    learner <- mlr3::lrn(
      "classif.ranger",
      predict_type = "prob",
      num.trees = as.integer(
        parameters$num.trees
      ),
      mtry.ratio = as.numeric(
        parameters$mtry.ratio
      ),
      min.node.size = as.integer(
        parameters$min.node.size
      ),
      num.threads = 1L,
      importance = "none",
      save.memory = TRUE,
      seed = as.integer(seed)
    )

    if (!is.null(class_weights)) {
      learner$param_set$values$class.weights <- class_weights
    }

    return(
      learner
    )
  }

  if (learner_id == "LightGBM") {
    return(
      mlr3::lrn(
        "classif.lightgbm",
        predict_type = "prob",
        learning_rate = as.numeric(
          parameters$learning_rate
        ),
        num_leaves = as.integer(
          parameters$num_leaves
        ),
        feature_fraction = as.numeric(
          parameters$feature_fraction
        ),
        bagging_fraction = as.numeric(
          parameters$bagging_fraction
        ),
        bagging_freq = 1L,
        min_data_in_leaf = as.integer(
          parameters$min_data_in_leaf
        ),
        num_iterations = as.integer(
          parameters$num_iterations
        ),
        num_threads = 1L,
        device_type = "cpu",
        deterministic = TRUE,
        seed = as.integer(seed)
      )
    )
  }

  if (learner_id == "XGBoost") {
    return(
      mlr3::lrn(
        "classif.xgboost",
        predict_type = "prob",
        eta = as.numeric(
          parameters$eta
        ),
        max_depth = as.integer(
          parameters$max_depth
        ),
        subsample = as.numeric(
          parameters$subsample
        ),
        colsample_bytree = as.numeric(
          parameters$colsample_bytree
        ),
        nrounds = as.integer(
          parameters$nrounds
        ),
        nthread = 1L,
        eval_metric = "logloss",
        seed = as.integer(seed)
      )
    )
  }

  if (learner_id == "RBF_SVM") {
    learner <- mlr3::lrn(
      "classif.svm",
      predict_type = "prob",
      type = "C-classification",
      kernel = "radial",
      cost = as.numeric(
        parameters$cost
      ),
      gamma = as.numeric(
        parameters$gamma
      ),
      scale = TRUE,
      cachesize = 500,
      tolerance = 1e-3,
      shrinking = TRUE
    )

    if (!is.null(class_weights)) {
      learner$param_set$values$class.weights <- class_weights
    }

    return(
      learner
    )
  }

  stop(
    "Unknown tunable learner: ",
    learner_id,
    call. = FALSE
  )
}

fit_predict_tunable_learner <- function(learner_id,
                                        parameters,
                                        train_matrix,
                                        train_metadata,
                                        validation_matrix,
                                        validation_metadata,
                                        seed,
                                        class_weights = NULL) {
  training_task <- make_mlr3_task(
    expression_matrix = train_matrix,
    metadata = train_metadata,
    task_id = paste0(
      "tune_train_",
      learner_id
    )
  )

  validation_task <- make_mlr3_task(
    expression_matrix = validation_matrix,
    metadata = validation_metadata,
    task_id = paste0(
      "tune_valid_",
      learner_id
    )
  )

  learner <- make_tunable_learner(
    learner_id = learner_id,
    parameters = parameters,
    seed = seed,
    class_weights = class_weights
  )

  set.seed(
    seed
  )

  learner$train(
    training_task$task
  )

  prediction <- learner$predict(
    validation_task$task
  )

  as.numeric(
    prediction$prob[
      ,
      POSITIVE_CLASS
    ]
  )
}


# ==============================================================================
# 10. MBO OPTIMISER
# ==============================================================================

make_mbo_optimizer <- function() {
  opt(
    "mbo",
    surrogate = srlrn(
      default_rf()
    )
  )
}

extract_archive_best_parameters <- function(instance,
                                            domain) {
  best <- instance$archive$best()

  if (
    "x_domain" %in% names(best) &&
      length(best$x_domain) > 0L
  ) {
    return(
      best$x_domain[[1L]]
    )
  }

  as.list(
    best[
      1L,
      domain$ids(),
      with = FALSE
    ]
  )
}


# ==============================================================================
# 11. CROSS-DISEASE RAW-DATA POOLING
# ==============================================================================

prepare_cross_disease_development <- function(disease_ids,
                                              common_genes) {
  component_counts <- list()
  component_metadata <- list()

  for (disease_id in disease_ids) {
    canonical <- load_canonical_disease(
      disease_id
    )

    counts <- canonical$development_counts[
      common_genes,
      ,
      drop = FALSE
    ]

    metadata <- canonical$development_metadata
    metadata$original_sample_id <- metadata$sample_id
    metadata$sample_id <- paste(
      disease_id,
      metadata$original_sample_id,
      sep = "__"
    )
    metadata$disease <- disease_id
    metadata$cohort_group <- factor(
      paste(
        disease_id,
        as.character(
          metadata$cohort_group
        ),
        sep = "__"
      )
    )

    colnames(
      counts
    ) <- metadata$sample_id

    component_counts[[disease_id]] <- counts
    component_metadata[[disease_id]] <- metadata
  }

  counts <- do.call(
    cbind,
    component_counts
  )

  metadata <- do.call(
    rbind,
    component_metadata
  )

  metadata <- as.data.frame(
    metadata,
    check.names = FALSE
  )

  metadata <- metadata[
    match(
      colnames(counts),
      metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  list(
    counts = counts,
    metadata = metadata
  )
}

load_cross_disease_common_genes <- function() {
  stage1_root <- dirname(
    get_config_value(
      config,
      "paths.stage_results.representation_specification",
      "results/01_representation_specification/TB"
    )
  )

  path <- file.path(
    stage1_root,
    "03_experiment_2_disease_composition",
    "cross_disease_common_genes.csv"
  )

  check_file_exists(
    path,
    "cross-disease common-gene manifest"
  )

  dt <- as.data.table(
    safe_read_csv(
      path
    )
  )

  check_required_columns(
    dt,
    "gene",
    "cross-disease common-gene manifest"
  )

  genes <- as.character(
    dt$gene
  )

  if (
    STRICT_AUDITS &&
      length(genes) != 4450L
  ) {
    stop(
      "Cross-disease common feature space contains ",
      length(genes),
      " genes; the dissertation design specifies 4,450.",
      call. = FALSE
    )
  }

  genes
}


# ==============================================================================
# 12. EXPERIMENT 1 - EQUAL-SIZE DEVELOPMENT BENCHMARK
# ==============================================================================

EXPERIMENT1_REPRESENTATIONS <- c(
  "VST_fSVA",
  "VST_ComBatRef",
  "ComBatSeq_log1p"
)

EXPERIMENT1_SELECTORS <- c(
  "information_gain",
  "importance",
  "boruta",
  "mrmr",
  "cmim"
)

EXPERIMENT1_LEARNERS <- c(
  "Featureless",
  "Naive_Bayes",
  "Linear_SVM",
  "RBF_SVM",
  "C5_0",
  "Random_Forest",
  "LightGBM",
  "XGBoost"
)

run_experiment_1 <- function() {
  output_root <- file.path(
    STAGE3_ROOT,
    "01_experiment_1_sample_support"
  )

  make_dir(
    output_root
  )

  experiment_diseases <- intersect(
    c(
      "TB",
      "RA",
      "UC",
      "CD"
    ),
    REQUESTED_DISEASES
  )

  for (disease_id in experiment_diseases) {
    message_section(
      paste0(
        "Stage 03 | Experiment 1 | ",
        disease_id
      )
    )

    data <- load_equal_size_development(
      disease_id
    )

    fold_map <- make_repeated_fold_map(
      metadata = data$metadata,
      folds = 10L,
      repeats = 2L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "experiment1"
      )
    )

    disease_dir <- file.path(
      output_root,
      disease_id
    )

    make_dir(
      disease_dir
    )

    safe_write_csv(
      fold_map,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_10fold_2repeat_fold_map.csv"
        )
      )
    )

    prediction_rows <- list()
    metric_rows <- list()
    gene_rows <- list()
    failure_rows <- list()

    prediction_counter <- 0L
    metric_counter <- 0L
    gene_counter <- 0L
    failure_counter <- 0L

    for (repeat_index in 1:2) {
      for (fold_id in 1:10) {
        validation_ids <- fold_map[`repeat` == repeat_index & fold == fold_id, sample_id]

        training_ids <- setdiff(
          data$metadata$sample_id,
          validation_ids
        )

        train_metadata <- data$metadata[
          match(
            training_ids,
            data$metadata$sample_id
          ),
          ,
          drop = FALSE
        ]

        validation_metadata <- data$metadata[
          match(
            validation_ids,
            data$metadata$sample_id
          ),
          ,
          drop = FALSE
        ]

        train_counts <- data$counts[
          ,
          train_metadata$sample_id,
          drop = FALSE
        ]

        validation_counts <- data$counts[
          ,
          validation_metadata$sample_id,
          drop = FALSE
        ]

        for (representation_index in seq_along(EXPERIMENT1_REPRESENTATIONS)) {
          representation_id <- EXPERIMENT1_REPRESENTATIONS[[representation_index]]

          fit_result <- tryCatch(
            stage1$fit_apply_representation(
              representation_id = representation_id,
              train_counts = train_counts,
              new_counts = validation_counts,
              train_metadata = train_metadata,
              new_metadata = validation_metadata,
              seed = make_trace_subseed(
                GLOBAL_SEED,
                disease_id,
                "experiment1",
                repeat_index,
                fold_id,
                representation_id
              )
            ),
            error = function(e) e
          )

          if (inherits(fit_result, "error")) {
            failure_counter <- failure_counter + 1L

            failure_rows[[failure_counter]] <- data.table(
              disease = disease_id,
              repeat = repeat_index,
              fold = fold_id,
              representation = representation_id,
              stage = "representation",
              selector = NA_character_,
              learner = NA_character_,
              error = conditionMessage(
                fit_result
              )
            )

            next
          }

          for (selector_index in seq_along(EXPERIMENT1_SELECTORS)) {
            selector_id <- EXPERIMENT1_SELECTORS[[selector_index]]

            selected <- tryCatch(
              stage2$select_ranked_features(
                train_matrix = fit_result$train,
                train_metadata = train_metadata,
                filter_id = selector_id,
                top_n = 50L,
                seed = make_trace_subseed(
                  GLOBAL_SEED,
                  disease_id,
                  "experiment1",
                  repeat_index,
                  fold_id,
                  representation_id,
                  selector_id
                ),
                rf_trees = 300L
              ),
              error = function(e) e
            )

            if (inherits(selected, "error")) {
              failure_counter <- failure_counter + 1L

              failure_rows[[failure_counter]] <- data.table(
                disease = disease_id,
                repeat = repeat_index,
                fold = fold_id,
                representation = representation_id,
                stage = "feature_selection",
                selector = selector_id,
                learner = NA_character_,
                error = conditionMessage(
                  selected
                )
              )

              next
            }

            genes <- selected$genes

            train_matrix <- fit_result$train[
              genes,
              ,
              drop = FALSE
            ]

            validation_matrix <- fit_result$new[
              genes,
              ,
              drop = FALSE
            ]

            gene_counter <- gene_counter + 1L

            gene_rows[[gene_counter]] <- data.table(
              disease = disease_id,
              repeat = repeat_index,
              fold = fold_id,
              representation = representation_id,
              selector = selector_id,
              rank = seq_along(
                genes
              ),
              gene = genes
            )

            for (learner_index in seq_along(EXPERIMENT1_LEARNERS)) {
              learner_id <- EXPERIMENT1_LEARNERS[[learner_index]]

              current_seed <- make_trace_subseed(
                GLOBAL_SEED,
                disease_id,
                "experiment1",
                repeat_index,
                fold_id,
                representation_id,
                selector_id,
                learner_id
              )

              probability <- tryCatch(
                fit_predict_fixed_learner(
                  learner_id = learner_id,
                  train_matrix = train_matrix,
                  train_metadata = train_metadata,
                  validation_matrix = validation_matrix,
                  validation_metadata = validation_metadata,
                  seed = current_seed
                ),
                error = function(e) e
              )

              if (inherits(probability, "error")) {
                failure_counter <- failure_counter + 1L

                failure_rows[[failure_counter]] <- data.table(
                  disease = disease_id,
                  repeat = repeat_index,
                  fold = fold_id,
                  representation = representation_id,
                  stage = "learner",
                  selector = selector_id,
                  learner = learner_id,
                  error = conditionMessage(
                    probability
                  )
                )

                next
              }

              metrics <- classification_metrics(
                truth = validation_metadata$outcome,
                probability = probability,
                threshold = 0.50
              )

              prediction_counter <- prediction_counter + 1L

              prediction_rows[[prediction_counter]] <- data.table(
                disease = disease_id,
                repeat = repeat_index,
                fold = fold_id,
                representation = representation_id,
                selector = selector_id,
                learner = learner_id,
                sample_id = validation_metadata$sample_id,
                cohort = as.character(
                  validation_metadata$cohort_group
                ),
                truth = as.character(
                  validation_metadata$outcome
                ),
                probability = probability
              )

              metric_counter <- metric_counter + 1L

              metric_rows[[metric_counter]] <- cbind(
                data.table(
                  disease = disease_id,
                  repeat = repeat_index,
                  fold = fold_id,
                  representation = representation_id,
                  selector = selector_id,
                  learner = learner_id
                ),
                metrics
              )
            }
          }
        }
      }
    }

    predictions <- rbindlist(
      prediction_rows,
      fill = TRUE
    )

    fold_metrics <- rbindlist(
      metric_rows,
      fill = TRUE
    )

    selected_genes <- rbindlist(
      gene_rows,
      fill = TRUE
    )

    failures <- rbindlist(
      failure_rows,
      fill = TRUE
    )

    safe_write_csv(
      predictions,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_prediction_level.csv"
        )
      )
    )

    safe_write_csv(
      fold_metrics,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_fold_metrics.csv"
        )
      )
    )

    safe_write_csv(
      selected_genes,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_fold_internal_selected_genes.csv"
        )
      )
    )

    if (nrow(failures) > 0L) {
      safe_write_csv(
        failures,
        file.path(
          disease_dir,
          paste0(
            disease_id,
            "_development_failures.csv"
          )
        )
      )
    }

    summary <- fold_metrics[
      ,
      .(
        successful_folds = .N,
        mean_AUROC = mean(
          AUROC,
          na.rm = TRUE
        ),
        sd_AUROC = stats::sd(
          AUROC,
          na.rm = TRUE
        ),
        worst_AUROC = min(
          AUROC,
          na.rm = TRUE
        ),
        mean_AUPRC = mean(
          AUPRC,
          na.rm = TRUE
        ),
        mean_BACC = mean(
          BACC,
          na.rm = TRUE
        ),
        mean_Brier = mean(
          Brier,
          na.rm = TRUE
        )
      ),
      by = .(
        representation,
        selector,
        learner
      )
    ]

    setorder(
      summary,
      -mean_AUROC,
      sd_AUROC,
      -worst_AUROC
    )

    safe_write_csv(
      summary,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_pipeline_summary.csv"
        )
      )
    )

    representation_sensitivity <- summary[
      selector == "importance" &
        learner == "LightGBM"
    ]

    safe_write_csv(
      representation_sensitivity,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_RFimportance_LightGBM_representation_sensitivity.csv"
        )
      )
    )

    if (nrow(summary) > 0L) {
      locked_development_winner <- summary[1L]

      safe_write_csv(
        locked_development_winner,
        file.path(
          disease_dir,
          paste0(
            disease_id,
            "_development_best_pipeline.csv"
          )
        )
      )
    }
  }

  invisible(TRUE)
}


# ==============================================================================
# 13. EXPERIMENT 2 - CROSS-DISEASE NESTED CV
# ==============================================================================

CROSS_DISEASE_OUTER_FOLDS <- 5L
CROSS_DISEASE_INNER_FOLDS <- 10L
CROSS_DISEASE_TUNING_EVALS <- 20L
CROSS_DISEASE_SMOTE_K <- 5L

CROSS_DISEASE_SELECTORS <- c(
  "mrmr",
  "cmim"
)

CROSS_DISEASE_CLASS_HANDLING <- c(
  "None",
  "SMOTE"
)

CROSS_DISEASE_LEARNERS <- c(
  "Random_Forest",
  "LightGBM",
  "XGBoost",
  "RBF_SVM"
)

cross_disease_fold_pipeline <- function(train_counts,
                                        train_metadata,
                                        validation_counts,
                                        validation_metadata,
                                        selector_id,
                                        class_handling,
                                        learner_id,
                                        parameters,
                                        seed) {
  representation <- stage1$fit_apply_vst_fsva(
    train_counts = train_counts,
    new_counts = validation_counts,
    train_metadata = train_metadata,
    seed = make_trace_subseed(
      seed,
      "representation"
    ),
    include_cohort_in_vst_design = FALSE
  )

  selected <- stage2$select_ranked_features(
    train_matrix = representation$train,
    train_metadata = train_metadata,
    filter_id = selector_id,
    top_n = 50L,
    seed = make_trace_subseed(
      seed,
      "feature_selection"
    ),
    rf_trees = 300L
  )

  genes <- selected$genes

  training_matrix <- representation$train[
    genes,
    ,
    drop = FALSE
  ]

  validation_matrix <- representation$new[
    genes,
    ,
    drop = FALSE
  ]

  balancing_audit <- NULL

  if (class_handling == "SMOTE") {
    smote_result <- stratified_smote(
      expression_matrix = training_matrix,
      metadata = train_metadata,
      k = CROSS_DISEASE_SMOTE_K,
      seed = make_trace_subseed(
        seed,
        "smote"
      ),
      include_disease = TRUE
    )

    training_matrix <- smote_result$matrix
    train_metadata <- smote_result$metadata
    balancing_audit <- smote_result$audit
  }

  probability <- fit_predict_tunable_learner(
    learner_id = learner_id,
    parameters = parameters,
    train_matrix = training_matrix,
    train_metadata = train_metadata,
    validation_matrix = validation_matrix,
    validation_metadata = validation_metadata,
    seed = make_trace_subseed(
      seed,
      "learner"
    )
  )

  list(
    probability = probability,
    genes = genes,
    balancing_audit = balancing_audit
  )
}

evaluate_cross_disease_parameters <- function(model_counts,
                                              model_metadata,
                                              outer_training_ids,
                                              selector_id,
                                              class_handling,
                                              learner_id,
                                              parameters,
                                              inner_fold_id,
                                              seed) {
  fold_rows <- list()

  for (inner_fold in seq_len(CROSS_DISEASE_INNER_FOLDS)) {
    inner_validation_ids <- outer_training_ids[
      inner_fold_id ==
        inner_fold
    ]

    inner_training_ids <- setdiff(
      outer_training_ids,
      inner_validation_ids
    )

    inner_train_metadata <- model_metadata[
      match(
        inner_training_ids,
        model_metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    inner_validation_metadata <- model_metadata[
      match(
        inner_validation_ids,
        model_metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    inner_train_counts <- model_counts[
      ,
      inner_train_metadata$sample_id,
      drop = FALSE
    ]

    inner_validation_counts <- model_counts[
      ,
      inner_validation_metadata$sample_id,
      drop = FALSE
    ]

    result <- tryCatch(
      cross_disease_fold_pipeline(
        train_counts = inner_train_counts,
        train_metadata = inner_train_metadata,
        validation_counts = inner_validation_counts,
        validation_metadata = inner_validation_metadata,
        selector_id = selector_id,
        class_handling = class_handling,
        learner_id = learner_id,
        parameters = parameters,
        seed = make_trace_subseed(
          seed,
          "inner_fold",
          inner_fold
        )
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      fold_rows[[inner_fold]] <- data.table(
        inner_fold = inner_fold,
        BACC = 0,
        AUROC = NA_real_,
        error = conditionMessage(
          result
        )
      )

      next
    }

    metrics <- classification_metrics(
      truth = inner_validation_metadata$outcome,
      probability = result$probability,
      threshold = 0.50
    )

    fold_rows[[inner_fold]] <- cbind(
      data.table(
        inner_fold = inner_fold,
        error = NA_character_
      ),
      metrics[
        ,
        .(
          BACC,
          AUROC
        )
      ]
    )
  }

  fold_dt <- rbindlist(
    fold_rows,
    fill = TRUE
  )

  mean_bacc <- if (
    nrow(fold_dt) ==
      CROSS_DISEASE_INNER_FOLDS &&
      all(
        is.finite(
          fold_dt$BACC
        )
      )
  ) {
    mean(
      fold_dt$BACC
    )
  } else {
    0
  }

  list(
    mean_bacc = mean_bacc,
    folds = fold_dt
  )
}

tune_cross_disease_pipeline <- function(model_counts,
                                        model_metadata,
                                        outer_training_ids,
                                        selector_id,
                                        class_handling,
                                        learner_id,
                                        inner_fold_id,
                                        context_id,
                                        seed,
                                        output_dir,
                                        evaluations = CROSS_DISEASE_TUNING_EVALS) {
  domain <- CROSS_DISEASE_SEARCH_SPACES[[learner_id]]

  history <- list()
  evaluation_id <- 0L

  objective_function <- function(xs) {
    evaluation_id <<- evaluation_id + 1L

    result <- evaluate_cross_disease_parameters(
      model_counts = model_counts,
      model_metadata = model_metadata,
      outer_training_ids = outer_training_ids,
      selector_id = selector_id,
      class_handling = class_handling,
      learner_id = learner_id,
      parameters = xs,
      inner_fold_id = inner_fold_id,
      seed = make_trace_subseed(
        seed,
        "evaluation",
        evaluation_id
      )
    )

    current_history <- copy(
      result$folds
    )

    current_history[
      ,
      `:=`(
        evaluation_id = evaluation_id,
        mean_inner_BACC = result$mean_bacc,
        context = context_id,
        selector = selector_id,
        class_handling = class_handling,
        learner = learner_id
      )
    ]

    for (parameter_name in names(xs)) {
      current_history[
        ,
        (parameter_name) := xs[[parameter_name]]
      ]
    }

    history[[evaluation_id]] <<- current_history

    list(
      loss = 1 -
        result$mean_bacc
    )
  }

  objective <- ObjectiveRFun$new(
    fun = objective_function,
    domain = domain,
    codomain = paradox::ps(
      loss = paradox::p_dbl(
        tags = "minimize"
      )
    )
  )

  instance <- OptimInstanceBatchSingleCrit$new(
    objective = objective,
    search_space = domain,
    terminator = trm(
      "evals",
      n_evals = as.integer(
        evaluations
      )
    )
  )

  set.seed(
    seed
  )

  optimiser <- make_mbo_optimizer()

  optimisation_error <- tryCatch(
    {
      optimiser$optimize(
        instance
      )
      NULL
    },
    error = function(e) conditionMessage(
      e
    )
  )

  if (!is.null(optimisation_error)) {
    stop(
      "MBO failed for ",
      context_id,
      ": ",
      optimisation_error,
      call. = FALSE
    )
  }

  make_dir(
    output_dir
  )

  archive_path <- file.path(
    output_dir,
    paste0(
      safe_file_name(
        context_id
      ),
      "_MBO_archive.rds"
    )
  )

  safe_save_rds(
    instance$archive$data,
    archive_path
  )

  archive_csv <- copy(
    instance$archive$data
  )

  list_columns <- names(
    archive_csv
  )[
    vapply(
      archive_csv,
      is.list,
      logical(1)
    )
  ]

  if (length(list_columns) > 0L) {
    archive_csv[
      ,
      (list_columns) := NULL
    ]
  }

  safe_write_csv(
    archive_csv,
    file.path(
      output_dir,
      paste0(
        safe_file_name(
          context_id
        ),
        "_MBO_archive.csv"
      )
    )
  )

  history_dt <- rbindlist(
    history,
    fill = TRUE
  )

  safe_write_csv(
    history_dt,
    file.path(
      output_dir,
      paste0(
        safe_file_name(
          context_id
        ),
        "_inner_history.csv"
      )
    )
  )

  best_parameters <- extract_archive_best_parameters(
    instance,
    domain
  )

  best <- instance$archive$best()

  list(
    parameters = best_parameters,
    best_inner_BACC = 1 -
      best$loss[[1L]],
    archive = instance$archive$data,
    history = history_dt
  )
}

fit_cross_disease_outer_fold <- function(model_counts,
                                         model_metadata,
                                         outer_training_ids,
                                         outer_validation_ids,
                                         selector_id,
                                         class_handling,
                                         learner_id,
                                         parameters,
                                         seed) {
  train_metadata <- model_metadata[
    match(
      outer_training_ids,
      model_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  validation_metadata <- model_metadata[
    match(
      outer_validation_ids,
      model_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  train_counts <- model_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]

  validation_counts <- model_counts[
    ,
    validation_metadata$sample_id,
    drop = FALSE
  ]

  result <- cross_disease_fold_pipeline(
    train_counts = train_counts,
    train_metadata = train_metadata,
    validation_counts = validation_counts,
    validation_metadata = validation_metadata,
    selector_id = selector_id,
    class_handling = class_handling,
    learner_id = learner_id,
    parameters = parameters,
    seed = seed
  )

  metrics <- classification_metrics(
    truth = validation_metadata$outcome,
    probability = result$probability,
    threshold = 0.50
  )

  list(
    probability = result$probability,
    metrics = metrics,
    genes = result$genes,
    validation_metadata = validation_metadata
  )
}

youden_threshold <- function(truth,
                             probability) {
  truth <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  roc_object <- pROC::roc(
    response = truth,
    predictor = probability,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  coords <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = c(
      "threshold",
      "sensitivity",
      "specificity"
    ),
    transpose = FALSE
  )

  threshold <- as.numeric(
    coords$threshold[[1L]]
  )

  threshold
}

run_cross_disease_model <- function(model_name,
                                    disease_ids,
                                    common_genes,
                                    output_root) {
  message_section(
    paste0(
      "Experiment 2 nested CV | ",
      model_name
    )
  )

  pooled <- prepare_cross_disease_development(
    disease_ids = disease_ids,
    common_genes = common_genes
  )

  model_counts <- pooled$counts
  model_metadata <- pooled$metadata

  outer_fold_id <- make_joint_stratified_folds(
    metadata = model_metadata,
    number_of_folds = CROSS_DISEASE_OUTER_FOLDS,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "experiment2",
      model_name,
      "outer"
    ),
    include_disease = TRUE
  )

  outer_map <- data.table(
    sample_id = model_metadata$sample_id,
    disease = model_metadata$disease,
    cohort = as.character(
      model_metadata$cohort_group
    ),
    outcome = as.character(
      model_metadata$outcome
    ),
    outer_fold = outer_fold_id
  )

  model_dir <- file.path(
    output_root,
    model_name
  )

  make_dir(
    model_dir
  )

  safe_write_csv(
    outer_map,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_outer_fold_map.csv"
      )
    )
  )

  outer_predictions <- list()
  outer_metrics <- list()
  best_parameter_rows <- list()
  selected_gene_rows <- list()

  prediction_counter <- 0L
  metric_counter <- 0L
  parameter_counter <- 0L
  gene_counter <- 0L

  pipeline_grid <- CJ(
    selector = CROSS_DISEASE_SELECTORS,
    class_handling = CROSS_DISEASE_CLASS_HANDLING,
    learner = CROSS_DISEASE_LEARNERS,
    unique = TRUE
  )

  safe_write_csv(
    pipeline_grid,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_pipeline_grid.csv"
      )
    )
  )

  # Use an explicit base-R loop to avoid data.table name shadowing.
  for (outer_fold_index in seq_len(CROSS_DISEASE_OUTER_FOLDS)) {
    outer_validation_ids <- outer_map$sample_id[
      outer_map$outer_fold ==
        outer_fold_index
    ]

    outer_training_ids <- setdiff(
      model_metadata$sample_id,
      outer_validation_ids
    )

    outer_training_metadata <- model_metadata[
      match(
        outer_training_ids,
        model_metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    inner_fold_id <- make_joint_stratified_folds(
      metadata = outer_training_metadata,
      number_of_folds = CROSS_DISEASE_INNER_FOLDS,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "experiment2",
        model_name,
        "outer",
        outer_fold_index,
        "inner"
      ),
      include_disease = TRUE
    )

    for (pipeline_index in seq_len(nrow(pipeline_grid))) {
      selector_id <- pipeline_grid$selector[[pipeline_index]]
      class_handling <- pipeline_grid$class_handling[[pipeline_index]]
      learner_id <- pipeline_grid$learner[[pipeline_index]]

      context_id <- paste(
        model_name,
        paste0(
          "outer",
          outer_fold_index
        ),
        selector_id,
        class_handling,
        learner_id,
        sep = "__"
      )

      tuning <- tune_cross_disease_pipeline(
        model_counts = model_counts,
        model_metadata = model_metadata,
        outer_training_ids = outer_training_ids,
        selector_id = selector_id,
        class_handling = class_handling,
        learner_id = learner_id,
        inner_fold_id = inner_fold_id,
        context_id = context_id,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          "experiment2",
          context_id
        ),
        output_dir = file.path(
          model_dir,
          "inner_tuning",
          paste0(
            "outer_",
            outer_fold_index
          )
        )
      )

      outer_result <- fit_cross_disease_outer_fold(
        model_counts = model_counts,
        model_metadata = model_metadata,
        outer_training_ids = outer_training_ids,
        outer_validation_ids = outer_validation_ids,
        selector_id = selector_id,
        class_handling = class_handling,
        learner_id = learner_id,
        parameters = tuning$parameters,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          "experiment2",
          context_id,
          "outer_fit"
        )
      )

      prediction_counter <- prediction_counter + 1L

      outer_predictions[[prediction_counter]] <- data.table(
        model = model_name,
        outer_fold = outer_fold_index,
        selector = selector_id,
        class_handling = class_handling,
        learner = learner_id,
        sample_id = outer_result$validation_metadata$sample_id,
        disease = outer_result$validation_metadata$disease,
        cohort = as.character(
          outer_result$validation_metadata$cohort_group
        ),
        truth = as.character(
          outer_result$validation_metadata$outcome
        ),
        probability = outer_result$probability
      )

      metric_counter <- metric_counter + 1L

      outer_metrics[[metric_counter]] <- cbind(
        data.table(
          model = model_name,
          outer_fold = outer_fold_index,
          selector = selector_id,
          class_handling = class_handling,
          learner = learner_id,
          best_inner_BACC = tuning$best_inner_BACC
        ),
        outer_result$metrics
      )

      parameter_counter <- parameter_counter + 1L

      parameter_row <- data.table(
        model = model_name,
        outer_fold = outer_fold_index,
        selector = selector_id,
        class_handling = class_handling,
        learner = learner_id,
        best_inner_BACC = tuning$best_inner_BACC
      )

      for (parameter_name in names(tuning$parameters)) {
        parameter_row[
          ,
          (parameter_name) := tuning$parameters[[parameter_name]]
        ]
      }

      best_parameter_rows[[parameter_counter]] <- parameter_row

      gene_counter <- gene_counter + 1L

      selected_gene_rows[[gene_counter]] <- data.table(
        model = model_name,
        outer_fold = outer_fold_index,
        selector = selector_id,
        class_handling = class_handling,
        learner = learner_id,
        rank = seq_along(
          outer_result$genes
        ),
        gene = outer_result$genes
      )
    }
  }

  predictions <- rbindlist(
    outer_predictions,
    fill = TRUE
  )

  metrics <- rbindlist(
    outer_metrics,
    fill = TRUE
  )

  best_parameters <- rbindlist(
    best_parameter_rows,
    fill = TRUE
  )

  selected_genes <- rbindlist(
    selected_gene_rows,
    fill = TRUE
  )

  safe_write_csv(
    predictions,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_nested_outer_predictions.csv"
      )
    )
  )

  safe_write_csv(
    metrics,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_nested_outer_metrics.csv"
      )
    )
  )

  safe_write_csv(
    best_parameters,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_outer_best_hyperparameters.csv"
      )
    )
  )

  safe_write_csv(
    selected_genes,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_outer_selected_genes.csv"
      )
    )
  )

  summary <- metrics[
    ,
    .(
      outer_folds = .N,
      mean_outer_AUROC = mean(
        AUROC,
        na.rm = TRUE
      ),
      sd_outer_AUROC = stats::sd(
        AUROC,
        na.rm = TRUE
      ),
      worst_outer_AUROC = min(
        AUROC,
        na.rm = TRUE
      ),
      mean_outer_AUPRC = mean(
        AUPRC,
        na.rm = TRUE
      ),
      mean_outer_BACC = mean(
        BACC,
        na.rm = TRUE
      ),
      mean_outer_Brier = mean(
        Brier,
        na.rm = TRUE
      ),
      mean_best_inner_BACC = mean(
        best_inner_BACC,
        na.rm = TRUE
      )
    ),
    by = .(
      selector,
      class_handling,
      learner
    )
  ]

  setorder(
    summary,
    -mean_outer_AUROC,
    sd_outer_AUROC,
    -worst_outer_AUROC
  )

  safe_write_csv(
    summary,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_nested_pipeline_summary.csv"
      )
    )
  )

  locked_pipeline <- list(
    selector = "mrmr",
    class_handling = "None",
    learner = "LightGBM"
  )

  locked_outer_predictions <- predictions[
    selector ==
      locked_pipeline$selector &
      class_handling ==
        locked_pipeline$class_handling &
      learner ==
        locked_pipeline$learner
  ]

  if (nrow(locked_outer_predictions) == 0L) {
    stop(
      model_name,
      " locked cross-disease pipeline has no outer predictions.",
      call. = FALSE
    )
  }

  threshold <- youden_threshold(
    truth = locked_outer_predictions$truth,
    probability = locked_outer_predictions$probability
  )

  full_development_inner_fold <- make_joint_stratified_folds(
    metadata = model_metadata,
    number_of_folds = CROSS_DISEASE_INNER_FOLDS,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "experiment2",
      model_name,
      "final_full_development_tuning"
    ),
    include_disease = TRUE
  )

  final_tuning <- tune_cross_disease_pipeline(
    model_counts = model_counts,
    model_metadata = model_metadata,
    outer_training_ids = model_metadata$sample_id,
    selector_id = locked_pipeline$selector,
    class_handling = locked_pipeline$class_handling,
    learner_id = locked_pipeline$learner,
    inner_fold_id = full_development_inner_fold,
    context_id = paste0(
      model_name,
      "__final_locked_mRMR_None_LightGBM"
    ),
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "experiment2",
      model_name,
      "final_locked_tuning"
    ),
    output_dir = file.path(
      model_dir,
      "final_locked_tuning"
    )
  )

  lock <- data.table(
    model = model_name,
    representation = "VST_fSVA",
    selector = locked_pipeline$selector,
    top_n = 50L,
    class_handling = locked_pipeline$class_handling,
    learner = locked_pipeline$learner,
    development_threshold = threshold,
    tuning_evaluations = CROSS_DISEASE_TUNING_EVALS
  )

  for (parameter_name in names(final_tuning$parameters)) {
    lock[
      ,
      (parameter_name) := final_tuning$parameters[[parameter_name]]
    ]
  }

  safe_write_csv(
    lock,
    file.path(
      model_dir,
      paste0(
        model_name,
        "_LOCKED_PIPELINE.csv"
      )
    )
  )

  list(
    summary = summary,
    predictions = predictions,
    lock = lock
  )
}

run_experiment_2 <- function() {
  if (!all(
    c(
      "TB",
      "RA",
      "UC",
      "CD",
      "T1D"
    ) %in% REQUESTED_DISEASES
  )) {
    stop(
      "Experiment 2 requires TB, RA, UC, CD and T1D.",
      call. = FALSE
    )
  }

  output_root <- file.path(
    STAGE3_ROOT,
    "02_experiment_2_disease_composition"
  )

  make_dir(
    output_root
  )

  common_genes <- load_cross_disease_common_genes()

  model_a <- run_cross_disease_model(
    model_name = "Model_A",
    disease_ids = c(
      "RA",
      "UC",
      "CD",
      "T1D"
    ),
    common_genes = common_genes,
    output_root = output_root
  )

  model_b <- run_cross_disease_model(
    model_name = "Model_B",
    disease_ids = c(
      "TB",
      "RA",
      "UC",
      "CD",
      "T1D"
    ),
    common_genes = common_genes,
    output_root = output_root
  )

  matched_summary <- merge(
    model_a$summary,
    model_b$summary,
    by = c(
      "selector",
      "class_handling",
      "learner"
    ),
    suffixes = c(
      "_ModelA",
      "_ModelB"
    ),
    all = TRUE
  )

  matched_summary[
    ,
    delta_mean_outer_AUROC :=
      mean_outer_AUROC_ModelB -
        mean_outer_AUROC_ModelA
  ]

  safe_write_csv(
    matched_summary,
    file.path(
      output_root,
      "Model_A_vs_Model_B_matched_pipeline_summary.csv"
    )
  )

  invisible(
    list(
      Model_A = model_a,
      Model_B = model_b
    )
  )
}


# ==============================================================================
# 14. EXPERIMENT 3 - DISEASE-SPECIFIC LOCO OPTIMISATION
# ==============================================================================

RBF_SCREEN_CONFIGS <- data.table(
  Config = c(
    "RBF_1",
    "RBF_2",
    "RBF_3",
    "RBF_4"
  ),
  cost = c(
    1,
    3.5345224,
    1.6746829,
    0.5606065
  ),
  gamma = c(
    0.0200,
    0.001012607,
    0.0003574437,
    0.007523222
  )
)

EXPERIMENT3_REPRESENTATIONS <- c(
  "VST_only",
  "log2_CPM_plus_1"
)

EXPERIMENT3_FEATURE_COUNTS <- c(
  25L,
  50L,
  100L
)

EXPERIMENT3_WEIGHTING <- c(
  "Unweighted",
  "Weighted"
)

EXPECTED_EXPERIMENT3_LOCKS <- data.table(
  disease = c(
    "TB",
    "RA",
    "UC",
    "CD",
    "T1D"
  ),
  representation = c(
    "log2_CPM_plus_1",
    "VST_only",
    "log2_CPM_plus_1",
    "log2_CPM_plus_1",
    "log2_CPM_plus_1"
  ),
  top_n = c(
    100L,
    25L,
    100L,
    100L,
    25L
  ),
  weighting = c(
    "Unweighted",
    "Unweighted",
    "Weighted",
    "Unweighted",
    "Weighted"
  )
)

class_weight_vector <- function(outcome) {
  outcome <- factor(
    as.character(outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  counts <- table(
    outcome
  )

  total <- sum(
    counts
  )

  weights <- total /
    (
      length(counts) *
        as.numeric(
          counts
        )
    )

  names(weights) <- names(
    counts
  )

  weights
}

fit_predict_disease_rf <- function(train_matrix,
                                   train_metadata,
                                   validation_matrix,
                                   weighting,
                                   seed) {
  x_train <- as.data.frame(
    t(
      train_matrix
    ),
    check.names = FALSE
  )

  x_validation <- as.data.frame(
    t(
      validation_matrix
    ),
    check.names = FALSE
  )

  safe_names <- make.unique(
    make.names(
      colnames(x_train),
      unique = TRUE
    )
  )

  colnames(x_train) <- safe_names
  colnames(x_validation) <- safe_names

  x_train$outcome <- factor(
    as.character(
      train_metadata$outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  class_weights <- if (
    weighting == "Weighted"
  ) {
    class_weight_vector(
      train_metadata$outcome
    )
  } else {
    NULL
  }

  set.seed(
    seed
  )

  model <- ranger::ranger(
    dependent.variable.name = "outcome",
    data = x_train,
    probability = TRUE,
    num.trees = 250L,
    mtry = max(
      1L,
      floor(
        sqrt(
          nrow(
            train_matrix
          )
        )
      )
    ),
    min.node.size = 5L,
    class.weights = class_weights,
    importance = "none",
    num.threads = 1L,
    seed = as.integer(seed)
  )

  prediction <- predict(
    model,
    data = x_validation
  )$predictions

  as.numeric(
    prediction[
      ,
      POSITIVE_CLASS
    ]
  )
}

fit_predict_disease_rbf <- function(train_matrix,
                                    train_metadata,
                                    validation_matrix,
                                    weighting,
                                    cost,
                                    gamma,
                                    seed) {
  x_train <- t(
    train_matrix
  )

  x_validation <- t(
    validation_matrix
  )

  outcome <- factor(
    as.character(
      train_metadata$outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  weights <- if (
    weighting == "Weighted"
  ) {
    class_weight_vector(
      outcome
    )
  } else {
    NULL
  }

  set.seed(
    seed
  )

  model <- e1071::svm(
    x = x_train,
    y = outcome,
    type = "C-classification",
    kernel = "radial",
    cost = cost,
    gamma = gamma,
    probability = TRUE,
    scale = TRUE,
    shrinking = TRUE,
    class.weights = weights
  )

  predicted <- predict(
    model,
    x_validation,
    probability = TRUE
  )

  probabilities <- attr(
    predicted,
    "probabilities"
  )

  as.numeric(
    probabilities[
      ,
      POSITIVE_CLASS
    ]
  )
}

find_experiment3_pair_directory <- function(disease_id,
                                            target_cohort,
                                            representation_id) {
  stage1_root <- dirname(
    get_config_value(
      config,
      "paths.stage_results.representation_specification",
      "results/01_representation_specification/TB"
    )
  )

  file.path(
    stage1_root,
    "04_experiment_3_cohort_composition",
    disease_id,
    "LOCO",
    target_cohort,
    representation_id
  )
}

load_stage1_pair_local <- function(directory) {
  development_path <- list.files(
    directory,
    pattern = "_development_matrix\\.rds$",
    full.names = TRUE
  )

  evaluation_path <- list.files(
    directory,
    pattern = "_evaluation_matrix\\.rds$",
    full.names = TRUE
  )

  development_metadata_path <- list.files(
    directory,
    pattern = "_development_metadata\\.csv$",
    full.names = TRUE
  )

  evaluation_metadata_path <- list.files(
    directory,
    pattern = "_evaluation_metadata\\.csv$",
    full.names = TRUE
  )

  if (
    length(development_path) != 1L ||
      length(evaluation_path) != 1L ||
      length(development_metadata_path) != 1L ||
      length(evaluation_metadata_path) != 1L
  ) {
    stop(
      "Expected one complete Stage-01 representation pair under ",
      directory,
      call. = FALSE
    )
  }

  list(
    development_matrix = coerce_numeric_matrix(
      readRDS(
        development_path[[1]]
      )
    ),
    evaluation_matrix = coerce_numeric_matrix(
      readRDS(
        evaluation_path[[1]]
      )
    ),
    development_metadata = as.data.frame(
      safe_read_csv(
        development_metadata_path[[1]]
      ),
      check.names = FALSE
    ),
    evaluation_metadata = as.data.frame(
      safe_read_csv(
        evaluation_metadata_path[[1]]
      ),
      check.names = FALSE
    )
  )
}

load_experiment3_genes <- function(disease_id,
                                   target_cohort,
                                   representation_id,
                                   top_n) {
  stage2_root <- dirname(
    get_config_value(
      config,
      "paths.stage_results.feature_set_characterisation",
      "results/02_feature_set_characterisation/TB"
    )
  )

  directory <- file.path(
    stage2_root,
    "04_experiment_3_cohort_composition",
    disease_id,
    "LOCO",
    target_cohort,
    representation_id
  )

  pattern <- paste0(
    "_RF_importance_top",
    top_n,
    "_genes\\.rds$"
  )

  files <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) != 1L) {
    stop(
      "Expected one Stage-02 gene vector for ",
      disease_id,
      " / ",
      target_cohort,
      " / ",
      representation_id,
      " / top ",
      top_n,
      ".",
      call. = FALSE
    )
  }

  readRDS(
    files[[1]]
  )
}

summarise_loco_candidate <- function(metrics,
                                     grouping_columns) {
  metrics[
    ,
    .(
      Mean_LOCO_AUROC = mean(
        AUROC,
        na.rm = TRUE
      ),
      Worst_cohort_AUROC = min(
        AUROC,
        na.rm = TRUE
      ),
      SD_LOCO_AUROC = stats::sd(
        AUROC,
        na.rm = TRUE
      ),
      Mean_LOCO_AUPRC = mean(
        AUPRC,
        na.rm = TRUE
      ),
      Mean_LOCO_BACC = mean(
        BACC,
        na.rm = TRUE
      ),
      Mean_LOCO_Brier = mean(
        Brier,
        na.rm = TRUE
      ),
      Cohorts = .N
    ),
    by = grouping_columns
  ]
}

select_loco_threshold <- function(predictions) {
  thresholds <- seq(
    0.05,
    0.95,
    by = 0.01
  )

  rows <- list()

  for (threshold_index in seq_along(thresholds)) {
    threshold <- thresholds[[threshold_index]]

    cohort_rows <- predictions[
      ,
      {
        current <- classification_metrics(
          truth = truth,
          probability = probability,
          threshold = threshold
        )

        .(
          BACC = current$BACC
        )
      },
      by = target_cohort
    ]

    rows[[threshold_index]] <- data.table(
      threshold = threshold,
      Mean_cohort_BACC = mean(
        cohort_rows$BACC,
        na.rm = TRUE
      ),
      Worst_cohort_BACC = min(
        cohort_rows$BACC,
        na.rm = TRUE
      ),
      Distance_from_0_50 = abs(
        threshold -
          0.50
      )
    )
  }

  table <- rbindlist(
    rows
  )

  setorder(
    table,
    -Mean_cohort_BACC,
    -Worst_cohort_BACC,
    Distance_from_0_50,
    threshold
  )

  list(
    threshold = table$threshold[[1L]],
    table = table
  )
}

run_experiment_3 <- function() {
  output_root <- file.path(
    STAGE3_ROOT,
    "03_experiment_3_cohort_composition"
  )

  make_dir(
    output_root
  )

  for (disease_id in REQUESTED_DISEASES) {
    message_section(
      paste0(
        "Stage 03 | Experiment 3 | ",
        disease_id
      )
    )

    canonical <- load_canonical_disease(
      disease_id
    )

    target_cohorts <- unique(
      as.character(
        canonical$development_metadata$cohort_group
      )
    )

    prediction_rows <- list()
    metric_rows <- list()

    prediction_counter <- 0L
    metric_counter <- 0L

    for (target_cohort in target_cohorts) {
      target_key <- safe_file_name(
        target_cohort
      )

      for (representation_id in EXPERIMENT3_REPRESENTATIONS) {
        pair_directory <- find_experiment3_pair_directory(
          disease_id = disease_id,
          target_cohort = target_key,
          representation_id = representation_id
        )

        pair <- load_stage1_pair_local(
          pair_directory
        )

        pair$development_metadata$outcome <- factor(
          as.character(
            pair$development_metadata$outcome
          ),
          levels = c(
            NEGATIVE_CLASS,
            POSITIVE_CLASS
          )
        )

        pair$evaluation_metadata$outcome <- factor(
          as.character(
            pair$evaluation_metadata$outcome
          ),
          levels = c(
            NEGATIVE_CLASS,
            POSITIVE_CLASS
          )
        )

        for (top_n in EXPERIMENT3_FEATURE_COUNTS) {
          genes <- load_experiment3_genes(
            disease_id = disease_id,
            target_cohort = target_key,
            representation_id = representation_id,
            top_n = top_n
          )

          train_matrix <- pair$development_matrix[
            genes,
            ,
            drop = FALSE
          ]

          validation_matrix <- pair$evaluation_matrix[
            genes,
            ,
            drop = FALSE
          ]

          for (weighting in EXPERIMENT3_WEIGHTING) {
            rf_probability <- fit_predict_disease_rf(
              train_matrix = train_matrix,
              train_metadata = pair$development_metadata,
              validation_matrix = validation_matrix,
              weighting = weighting,
              seed = make_trace_subseed(
                GLOBAL_SEED,
                disease_id,
                "experiment3",
                target_cohort,
                representation_id,
                top_n,
                weighting,
                "RF"
              )
            )

            rf_metrics <- classification_metrics(
              truth = pair$evaluation_metadata$outcome,
              probability = rf_probability,
              threshold = 0.50
            )

            prediction_counter <- prediction_counter + 1L

            prediction_rows[[prediction_counter]] <- data.table(
              disease = disease_id,
              representation = representation_id,
              target_cohort = target_cohort,
              weighting = weighting,
              top_n = top_n,
              model = "Random_Forest",
              config = "RF_1",
              sample_id = pair$evaluation_metadata$sample_id,
              truth = as.character(
                pair$evaluation_metadata$outcome
              ),
              probability = rf_probability
            )

            metric_counter <- metric_counter + 1L

            metric_rows[[metric_counter]] <- cbind(
              data.table(
                disease = disease_id,
                representation = representation_id,
                target_cohort = target_cohort,
                weighting = weighting,
                top_n = top_n,
                model = "Random_Forest",
                config = "RF_1",
                cost = NA_real_,
                gamma = NA_real_
              ),
              rf_metrics
            )

            for (config_index in seq_len(nrow(RBF_SCREEN_CONFIGS))) {
              rbf_row <- RBF_SCREEN_CONFIGS[
                config_index
              ]

              rbf_probability <- fit_predict_disease_rbf(
                train_matrix = train_matrix,
                train_metadata = pair$development_metadata,
                validation_matrix = validation_matrix,
                weighting = weighting,
                cost = rbf_row$cost,
                gamma = rbf_row$gamma,
                seed = make_trace_subseed(
                  GLOBAL_SEED,
                  disease_id,
                  "experiment3",
                  target_cohort,
                  representation_id,
                  top_n,
                  weighting,
                  rbf_row$Config
                )
              )

              rbf_metrics <- classification_metrics(
                truth = pair$evaluation_metadata$outcome,
                probability = rbf_probability,
                threshold = 0.50
              )

              prediction_counter <- prediction_counter + 1L

              prediction_rows[[prediction_counter]] <- data.table(
                disease = disease_id,
                representation = representation_id,
                target_cohort = target_cohort,
                weighting = weighting,
                top_n = top_n,
                model = "RBF_SVM",
                config = rbf_row$Config,
                sample_id = pair$evaluation_metadata$sample_id,
                truth = as.character(
                  pair$evaluation_metadata$outcome
                ),
                probability = rbf_probability
              )

              metric_counter <- metric_counter + 1L

              metric_rows[[metric_counter]] <- cbind(
                data.table(
                  disease = disease_id,
                  representation = representation_id,
                  target_cohort = target_cohort,
                  weighting = weighting,
                  top_n = top_n,
                  model = "RBF_SVM",
                  config = rbf_row$Config,
                  cost = rbf_row$cost,
                  gamma = rbf_row$gamma
                ),
                rbf_metrics
              )
            }
          }
        }
      }
    }

    predictions <- rbindlist(
      prediction_rows,
      fill = TRUE
    )

    metrics <- rbindlist(
      metric_rows,
      fill = TRUE
    )

    disease_dir <- file.path(
      output_root,
      disease_id
    )

    make_dir(
      disease_dir
    )

    safe_write_csv(
      predictions,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_LOCO_predictions_all_candidates.csv"
        )
      )
    )

    safe_write_csv(
      metrics,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_LOCO_metrics_all_candidates.csv"
        )
      )
    )

    screen_summary <- summarise_loco_candidate(
      metrics,
      c(
        "representation",
        "weighting",
        "top_n",
        "model",
        "config",
        "cost",
        "gamma"
      )
    )

    setorder(
      screen_summary,
      -Worst_cohort_AUROC,
      -Mean_LOCO_AUROC,
      SD_LOCO_AUROC,
      -Mean_LOCO_BACC
    )

    safe_write_csv(
      screen_summary,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LOCO_initial_screen_summary.csv"
        )
      )
    )

    learner_best <- screen_summary[
      ,
      .SD[1L],
      by = model
    ]

    setorder(
      learner_best,
      -Worst_cohort_AUROC,
      -Mean_LOCO_AUROC,
      SD_LOCO_AUROC,
      -Mean_LOCO_BACC
    )

    selected_learner <- learner_best$model[[1L]]

    if (
      selected_learner !=
        "Random_Forest"
    ) {
      warning(
        disease_id,
        " initial LOCO screen selected ",
        selected_learner,
        " by the robust hierarchy; the final dissertation reports Random Forest convergence. ",
        "Inspect the Stage-01/02 inputs before using the lock.",
        call. = FALSE
      )
    }

    rf_summary <- screen_summary[
      model ==
        "Random_Forest"
    ]

    representation_best <- rf_summary[
      ,
      .SD[
        order(
          -Worst_cohort_AUROC,
          -Mean_LOCO_AUROC,
          SD_LOCO_AUROC,
          -Mean_LOCO_BACC
        )
      ][1L],
      by = representation
    ]

    setorder(
      representation_best,
      -Worst_cohort_AUROC,
      -Mean_LOCO_AUROC,
      SD_LOCO_AUROC,
      -Mean_LOCO_BACC
    )

    selected_representation <- representation_best$representation[[1L]]

    definitive <- rf_summary[
      representation ==
        selected_representation
    ]

    setorder(
      definitive,
      -Worst_cohort_AUROC,
      -Mean_LOCO_AUROC,
      SD_LOCO_AUROC,
      -Mean_LOCO_BACC
    )

    locked <- definitive[1L]

    locked_predictions <- predictions[
      representation ==
        locked$representation &
      weighting ==
        locked$weighting &
      top_n ==
        locked$top_n &
      model ==
        "Random_Forest"
    ]

    threshold_result <- select_loco_threshold(
      locked_predictions
    )

    locked[
      ,
      development_threshold :=
        threshold_result$threshold
    ]

    safe_write_csv(
      definitive,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_six_RF_candidates_selected_representation.csv"
        )
      )
    )

    safe_write_csv(
      threshold_result$table,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_development_LOCO_threshold_search.csv"
        )
      )
    )

    safe_write_csv(
      locked,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LOCKED_DISEASE_SPECIFIC_PIPELINE.csv"
        )
      )
    )

    expected <- EXPECTED_EXPERIMENT3_LOCKS[
      disease ==
        disease_id
    ]

    if (nrow(expected) == 1L) {
      matches_expected <- (
        locked$representation ==
          expected$representation &&
          locked$top_n ==
            expected$top_n &&
          locked$weighting ==
            expected$weighting
      )

      audit <- data.table(
        disease = disease_id,
        observed_representation = locked$representation,
        expected_representation = expected$representation,
        observed_top_n = locked$top_n,
        expected_top_n = expected$top_n,
        observed_weighting = locked$weighting,
        expected_weighting = expected$weighting,
        matches_dissertation = matches_expected
      )

      safe_write_csv(
        audit,
        file.path(
          disease_dir,
          paste0(
            disease_id,
            "_lock_reproduction_audit.csv"
          )
        )
      )

      if (
        STRICT_AUDITS &&
          !matches_expected
      ) {
        stop(
          disease_id,
          " Stage-03 lock differs from the dissertation result.",
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}


# ==============================================================================
# 15. EXPERIMENT 4 - PYTHON INITIALISATION
# ==============================================================================

initialise_architecture_python <- function() {
  require_packages(
    "reticulate",
    bioc = FALSE
  )

  python_override <- Sys.getenv(
    "TRACE_OMICS_PYTHON",
    unset = ""
  )

  if (nzchar(
    python_override
  )) {
    reticulate::use_python(
      python_override,
      required = TRUE
    )
  }

  np <- reticulate::import(
    "numpy",
    delay_load = FALSE,
    convert = FALSE
  )

  lightgbm_py <- reticulate::import(
    "lightgbm",
    delay_load = FALSE,
    convert = FALSE
  )

  tabpfn_py <- reticulate::import(
    "tabpfn",
    delay_load = FALSE,
    convert = FALSE
  )

  torch_py <- reticulate::import(
    "torch",
    delay_load = FALSE,
    convert = FALSE
  )

  cuda_available <- tryCatch(
    isTRUE(
      reticulate::py_to_r(
        torch_py$cuda$is_available()
      )
    ),
    error = function(e) FALSE
  )

  device <- if (
    cuda_available
  ) {
    "cuda"
  } else {
    "cpu"
  }

  versions <- data.table(
    component = c(
      "Python",
      "NumPy",
      "LightGBM",
      "TabPFN",
      "PyTorch"
    ),
    observed = c(
      reticulate::py_config()$version_string,
      as.character(
        reticulate::py_to_r(
          np$`__version__`
        )
      ),
      as.character(
        reticulate::py_to_r(
          lightgbm_py$`__version__`
        )
      ),
      as.character(
        reticulate::py_to_r(
          tabpfn_py$`__version__`
        )
      ),
      as.character(
        reticulate::py_to_r(
          torch_py$`__version__`
        )
      )
    ),
    expected = c(
      "3.10",
      "1.22.4",
      "4.6.0",
      "0.1.11",
      "2.4.1"
    )
  )

  safe_write_csv(
    versions,
    file.path(
      AUDIT_ROOT,
      "Experiment4_python_version_audit.csv"
    )
  )

  list(
    np = np,
    lightgbm = lightgbm_py,
    tabpfn = tabpfn_py,
    torch = torch_py,
    device = device,
    versions = versions
  )
}


# ==============================================================================
# 16. EXPERIMENT 4 - LIGHTGBM GRID AND PYTHON FIT/PREDICT
# ==============================================================================

LIGHTGBM_GRID <- data.table(
  Config_ID = paste0(
    "LGBM_",
    sprintf(
      "%02d",
      1:10
    )
  ),
  n_estimators = c(
    100L,
    200L,
    200L,
    300L,
    300L,
    150L,
    250L,
    400L,
    250L,
    500L
  ),
  learning_rate = c(
    0.05,
    0.03,
    0.05,
    0.03,
    0.02,
    0.10,
    0.05,
    0.02,
    0.03,
    0.015
  ),
  num_leaves = c(
    7L,
    7L,
    15L,
    15L,
    15L,
    7L,
    31L,
    15L,
    31L,
    15L
  ),
  max_depth = c(
    3L,
    3L,
    4L,
    4L,
    5L,
    3L,
    5L,
    5L,
    6L,
    4L
  ),
  min_child_samples = c(
    10L,
    10L,
    10L,
    15L,
    20L,
    20L,
    15L,
    20L,
    20L,
    10L
  ),
  subsample = c(
    0.80,
    0.80,
    0.80,
    0.90,
    0.80,
    0.80,
    0.80,
    1.00,
    0.80,
    0.90
  ),
  colsample_bytree = c(
    0.80,
    0.80,
    0.80,
    0.90,
    0.80,
    0.80,
    0.80,
    1.00,
    0.90,
    0.80
  ),
  reg_alpha = c(
    0,
    0,
    0,
    0,
    0.5,
    0,
    0.5,
    0,
    1,
    0.5
  ),
  reg_lambda = c(
    1,
    1,
    1,
    2,
    2,
    1,
    2,
    5,
    5,
    2
  )
)

to_numpy_x <- function(x,
                       python) {
  python$np$array(
    x,
    dtype = "float32"
  )
}

to_numpy_y <- function(y,
                       python) {
  y_numeric <- as.integer(
    as.character(
      y
    ) ==
      POSITIVE_CLASS
  )

  python$np$array(
    y_numeric,
    dtype = "int64"
  )
}

fit_predict_python_lightgbm <- function(x_train,
                                        y_train,
                                        x_test,
                                        parameters,
                                        seed,
                                        python) {
  model <- python$lightgbm$LGBMClassifier(
    objective = "binary",
    n_estimators = as.integer(
      parameters$n_estimators[[1L]]
    ),
    learning_rate = parameters$learning_rate[[1L]],
    num_leaves = as.integer(
      parameters$num_leaves[[1L]]
    ),
    max_depth = as.integer(
      parameters$max_depth[[1L]]
    ),
    min_child_samples = as.integer(
      parameters$min_child_samples[[1L]]
    ),
    subsample = parameters$subsample[[1L]],
    subsample_freq = 1L,
    colsample_bytree = parameters$colsample_bytree[[1L]],
    reg_alpha = parameters$reg_alpha[[1L]],
    reg_lambda = parameters$reg_lambda[[1L]],
    class_weight = "balanced",
    random_state = as.integer(
      seed
    ),
    n_jobs = -1L,
    verbosity = -1L
  )

  model$fit(
    to_numpy_x(
      x_train,
      python
    ),
    to_numpy_y(
      y_train,
      python
    )
  )

  probability_matrix <- reticulate::py_to_r(
    model$predict_proba(
      to_numpy_x(
        x_test,
        python
      )
    )
  )

  as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )
}

fit_predict_python_tabpfn <- function(x_train,
                                      y_train,
                                      x_test,
                                      seed,
                                      python) {
  seed_integer <- as.integer(
    seed %%
      .Machine$integer.max
  )

  python$np$random$seed(
    seed_integer
  )

  model <- python$tabpfn$TabPFNClassifier(
    device = python$device,
    seed = seed_integer
  )

  model$fit(
    to_numpy_x(
      x_train,
      python
    ),
    to_numpy_y(
      y_train,
      python
    )
  )

  probability_matrix <- reticulate::py_to_r(
    model$predict_proba(
      to_numpy_x(
        x_test,
        python
      )
    )
  )

  as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )
}


# ==============================================================================
# 17. EXPERIMENT 4 - OUTER/INNER ARCHITECTURE DEVELOPMENT
# ==============================================================================

load_architecture_outer_fold_pair <- function(disease_id,
                                              fold_id) {
  stage2_root <- dirname(
    get_config_value(
      config,
      "paths.stage_results.feature_set_characterisation",
      "results/02_feature_set_characterisation/TB"
    )
  )

  directory <- file.path(
    stage2_root,
    "05_experiment_4_model_architecture",
    disease_id,
    "outer_cv",
    paste0(
      "fold_",
      fold_id
    )
  )

  train_path <- list.files(
    directory,
    pattern = "_development_selected_matrix\\.rds$",
    full.names = TRUE
  )

  validation_path <- list.files(
    directory,
    pattern = "_evaluation_selected_matrix\\.rds$",
    full.names = TRUE
  )

  genes_path <- list.files(
    directory,
    pattern = "_top50_genes\\.rds$",
    full.names = TRUE
  )

  if (
    length(train_path) != 1L ||
      length(validation_path) != 1L ||
      length(genes_path) != 1L
  ) {
    stop(
      "Incomplete Stage-02 architecture fold bundle: ",
      directory,
      call. = FALSE
    )
  }

  list(
    train = coerce_numeric_matrix(
      readRDS(
        train_path[[1]]
      )
    ),
    validation = coerce_numeric_matrix(
      readRDS(
        validation_path[[1]]
      )
    ),
    genes = readRDS(
      genes_path[[1]]
    )
  )
}

architecture_inner_tuning <- function(disease_id,
                                      outer_training_ids,
                                      canonical,
                                      python,
                                      seed,
                                      outer_fold) {
  outer_metadata <- canonical$development_metadata[
    match(
      outer_training_ids,
      canonical$development_metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  inner_fold_id <- make_joint_stratified_folds(
    metadata = outer_metadata,
    number_of_folds = 3L,
    seed = make_trace_subseed(
      seed,
      "architecture_inner"
    )
  )

  rows <- list()
  counter <- 0L

  for (inner_fold in 1:3) {
    validation_ids <- outer_training_ids[
      inner_fold_id ==
        inner_fold
    ]

    training_ids <- setdiff(
      outer_training_ids,
      validation_ids
    )

    train_metadata <- canonical$development_metadata[
      match(
        training_ids,
        canonical$development_metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    validation_metadata <- canonical$development_metadata[
      match(
        validation_ids,
        canonical$development_metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    train_counts <- canonical$development_counts[
      ,
      train_metadata$sample_id,
      drop = FALSE
    ]

    validation_counts <- canonical$development_counts[
      ,
      validation_metadata$sample_id,
      drop = FALSE
    ]

    representation <- stage2$fit_apply_log2_cpm_feature_stage(
      train_counts = train_counts,
      validation_counts = validation_counts,
      min_total_count = as.integer(
        get_config_value(
          config,
          "analysis.min_gene_count",
          10
        )
      ),
      min_features = 50L
    )

    selected <- stage2$select_ranked_features(
      train_matrix = representation$train,
      train_metadata = train_metadata,
      filter_id = "importance",
      top_n = 50L,
      seed = make_trace_subseed(
        seed,
        "inner",
        inner_fold,
        "RF_importance"
      ),
      rf_trees = 250L
    )

    genes <- selected$genes

    x_train <- t(
      representation$train[
        genes,
        ,
        drop = FALSE
      ]
    )

    x_validation <- t(
      representation$validation[
        genes,
        ,
        drop = FALSE
      ]
    )

    for (config_index in seq_len(nrow(LIGHTGBM_GRID))) {
      config_row <- LIGHTGBM_GRID[
        config_index
      ]

      probability <- fit_predict_python_lightgbm(
        x_train = x_train,
        y_train = train_metadata$outcome,
        x_test = x_validation,
        parameters = config_row,
        seed = make_trace_subseed(
          seed,
          "inner",
          inner_fold,
          config_row$Config_ID
        ),
        python = python
      )

      counter <- counter + 1L

      rows[[counter]] <- data.table(
        disease = disease_id,
        outer_fold = outer_fold,
        inner_fold = inner_fold,
        Config_ID = config_row$Config_ID,
        AUROC = safe_auroc(
          validation_metadata$outcome,
          probability
        )
      )
    }
  }

  raw <- rbindlist(
    rows
  )

  summary <- raw[
    ,
    .(
      Mean_inner_AUROC = mean(
        AUROC,
        na.rm = TRUE
      ),
      SD_inner_AUROC = stats::sd(
        AUROC,
        na.rm = TRUE
      )
    ),
    by = Config_ID
  ]

  setorder(
    summary,
    -Mean_inner_AUROC,
    SD_inner_AUROC,
    Config_ID
  )

  list(
    raw = raw,
    summary = summary,
    best_config_id = summary$Config_ID[[1L]]
  )
}

run_experiment_4 <- function() {
  output_root <- file.path(
    STAGE3_ROOT,
    "04_experiment_4_model_architecture"
  )

  make_dir(
    output_root
  )

  python <- initialise_architecture_python()

  safe_write_csv(
    LIGHTGBM_GRID,
    file.path(
      output_root,
      "prespecified_LightGBM_10_configuration_grid.csv"
    )
  )

  for (disease_id in REQUESTED_DISEASES) {
    message_section(
      paste0(
        "Stage 03 | Experiment 4 | ",
        disease_id
      )
    )

    canonical <- load_canonical_disease(
      disease_id
    )

    stage1_root <- dirname(
      get_config_value(
        config,
        "paths.stage_results.representation_specification",
        "results/01_representation_specification/TB"
      )
    )

    fold_map_path <- file.path(
      stage1_root,
      "05_experiment_4_model_architecture",
      disease_id,
      paste0(
        disease_id,
        "_architecture_outer_fold_map.csv"
      )
    )

    check_file_exists(
      fold_map_path,
      paste0(
        disease_id,
        " architecture outer fold map"
      )
    )

    fold_map <- as.data.table(
      safe_read_csv(
        fold_map_path
      )
    )

    prediction_rows <- list()
    metric_rows <- list()
    tuning_rows <- list()
    tuning_raw_rows <- list()
    gene_rows <- list()

    for (outer_fold in 1:5) {
      validation_ids <- fold_map[
        fold ==
          outer_fold,
        sample_id
      ]

      training_ids <- setdiff(
        canonical$development_metadata$sample_id,
        validation_ids
      )

      tuning <- architecture_inner_tuning(
        disease_id = disease_id,
        outer_training_ids = training_ids,
        canonical = canonical,
        python = python,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment4",
          "outer",
          outer_fold
        ),
        outer_fold = outer_fold
      )

      tuning_rows[[outer_fold]] <- tuning$summary

      tuning_rows[[outer_fold]][
        ,
        outer_fold := outer_fold
      ]

      tuning_raw_rows[[outer_fold]] <- tuning$raw

      outer_pair <- load_architecture_outer_fold_pair(
        disease_id = disease_id,
        fold_id = outer_fold
      )

      if (
        !identical(
          colnames(
            outer_pair$train
          ),
          training_ids
        )
      ) {
        if (!setequal(
          colnames(
            outer_pair$train
          ),
          training_ids
        )) {
          stop(
            disease_id,
            " architecture outer-training sample mismatch in fold ",
            outer_fold,
            ".",
            call. = FALSE
          )
        }

        outer_pair$train <- outer_pair$train[
          ,
          training_ids,
          drop = FALSE
        ]
      }

      if (
        !identical(
          colnames(
            outer_pair$validation
          ),
          validation_ids
        )
      ) {
        if (!setequal(
          colnames(
            outer_pair$validation
          ),
          validation_ids
        )) {
          stop(
            disease_id,
            " architecture outer-validation sample mismatch in fold ",
            outer_fold,
            ".",
            call. = FALSE
          )
        }

        outer_pair$validation <- outer_pair$validation[
          ,
          validation_ids,
          drop = FALSE
        ]
      }

      train_metadata <- canonical$development_metadata[
        match(
          training_ids,
          canonical$development_metadata$sample_id
        ),
        ,
        drop = FALSE
      ]

      validation_metadata <- canonical$development_metadata[
        match(
          validation_ids,
          canonical$development_metadata$sample_id
        ),
        ,
        drop = FALSE
      ]

      best_config <- LIGHTGBM_GRID[
        Config_ID ==
          tuning$best_config_id
      ]

      x_train <- t(
        outer_pair$train
      )

      x_validation <- t(
        outer_pair$validation
      )

      lightgbm_probability <- fit_predict_python_lightgbm(
        x_train = x_train,
        y_train = train_metadata$outcome,
        x_test = x_validation,
        parameters = best_config,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment4",
          outer_fold,
          "LightGBM"
        ),
        python = python
      )

      tabpfn_probability <- fit_predict_python_tabpfn(
        x_train = x_train,
        y_train = train_metadata$outcome,
        x_test = x_validation,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment4",
          outer_fold,
          "TabPFN"
        ),
        python = python
      )

      for (model_name in c(
        "LightGBM",
        "TabPFN"
      )) {
        probability <- if (
          model_name ==
            "LightGBM"
        ) {
          lightgbm_probability
        } else {
          tabpfn_probability
        }

        prediction_rows[[length(prediction_rows) + 1L]] <- data.table(
          disease = disease_id,
          outer_fold = outer_fold,
          model = model_name,
          selected_LightGBM_config = tuning$best_config_id,
          sample_id = validation_metadata$sample_id,
          cohort = as.character(
            validation_metadata$cohort_group
          ),
          truth = as.character(
            validation_metadata$outcome
          ),
          probability = probability
        )

        metrics <- classification_metrics(
          truth = validation_metadata$outcome,
          probability = probability,
          threshold = 0.50
        )

        metric_rows[[length(metric_rows) + 1L]] <- cbind(
          data.table(
            disease = disease_id,
            outer_fold = outer_fold,
            model = model_name,
            selected_LightGBM_config = tuning$best_config_id
          ),
          metrics
        )
      }

      gene_rows[[outer_fold]] <- data.table(
        disease = disease_id,
        outer_fold = outer_fold,
        rank = seq_along(
          outer_pair$genes
        ),
        gene = outer_pair$genes
      )
    }

    predictions <- rbindlist(
      prediction_rows,
      fill = TRUE
    )

    metrics <- rbindlist(
      metric_rows,
      fill = TRUE
    )

    tuning_summary_by_outer <- rbindlist(
      tuning_rows,
      fill = TRUE
    )

    tuning_raw <- rbindlist(
      tuning_raw_rows,
      fill = TRUE
    )

    genes <- rbindlist(
      gene_rows,
      fill = TRUE
    )

    disease_dir <- file.path(
      output_root,
      disease_id
    )

    make_dir(
      disease_dir
    )

    safe_write_csv(
      predictions,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_architecture_outer_predictions.csv"
        )
      )
    )

    safe_write_csv(
      metrics,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_architecture_outer_metrics.csv"
        )
      )
    )

    safe_write_csv(
      tuning_summary_by_outer,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LightGBM_inner_tuning_by_outer_fold.csv"
        )
      )
    )

    safe_write_csv(
      tuning_raw,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LightGBM_inner_tuning_raw.csv"
        )
      )
    )

    safe_write_csv(
      genes,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_architecture_outer_selected_genes.csv"
        )
      )
    )

    model_summary <- metrics[
      ,
      .(
        Outer_folds = .N,
        Mean_AUROC = mean(
          AUROC,
          na.rm = TRUE
        ),
        SD_AUROC = stats::sd(
          AUROC,
          na.rm = TRUE
        ),
        Worst_AUROC = min(
          AUROC,
          na.rm = TRUE
        ),
        Best_AUROC = max(
          AUROC,
          na.rm = TRUE
        ),
        Mean_AUPRC = mean(
          AUPRC,
          na.rm = TRUE
        ),
        Mean_BACC = mean(
          BACC,
          na.rm = TRUE
        ),
        Mean_Brier = mean(
          Brier,
          na.rm = TRUE
        )
      ),
      by = model
    ]

    safe_write_csv(
      model_summary,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_architecture_outer_summary.csv"
        )
      )
    )

    final_lgbm_tuning_summary <- tuning_summary_by_outer[
      ,
      .(
        Outer_folds = .N,
        Mean_inner_AUROC = mean(
          Mean_inner_AUROC,
          na.rm = TRUE
        ),
        SD_inner_AUROC = stats::sd(
          Mean_inner_AUROC,
          na.rm = TRUE
        ),
        Worst_outer_inner_AUROC = min(
          Mean_inner_AUROC,
          na.rm = TRUE
        )
      ),
      by = Config_ID
    ]

    setorder(
      final_lgbm_tuning_summary,
      -Mean_inner_AUROC,
      SD_inner_AUROC,
      -Worst_outer_inner_AUROC,
      Config_ID
    )

    final_config_id <- final_lgbm_tuning_summary$Config_ID[[1L]]

    final_config <- LIGHTGBM_GRID[
      Config_ID ==
        final_config_id
    ]

    safe_write_csv(
      final_lgbm_tuning_summary,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_final_LightGBM_config_selection.csv"
        )
      )
    )

    stage2_root <- dirname(
      get_config_value(
        config,
        "paths.stage_results.feature_set_characterisation",
        "results/02_feature_set_characterisation/TB"
      )
    )

    full_gene_directory <- file.path(
      stage2_root,
      "05_experiment_4_model_architecture",
      disease_id,
      "frozen_diagnostic_test"
    )

    final_gene_candidates <- list.files(
      full_gene_directory,
      pattern = "_architecture_top50_genes\\.rds$",
      full.names = TRUE
    )

    if (length(final_gene_candidates) != 1L) {
      stop(
        "Expected one development-derived architecture top-50 gene vector for ",
        disease_id,
        ".",
        call. = FALSE
      )
    }

    final_genes <- readRDS(
      final_gene_candidates[[1]]
    )

    if (length(final_genes) != 50L) {
      stop(
        disease_id,
        " architecture lock contains ",
        length(final_genes),
        " genes; expected 50.",
        call. = FALSE
      )
    }

    final_lock <- data.table(
      disease = disease_id,
      representation = "log2_CPM_plus_1",
      feature_selection = "RF permutation importance",
      top_n = 50L,
      LightGBM_config_id = final_config_id,
      TabPFN_version = "0.1.11",
      TabPFN_device = python$device,
      development_gene_vector_source = final_gene_candidates[[1]]
    )

    final_lock <- cbind(
      final_lock,
      final_config[
        ,
        !c(
          "Config_ID"
        ),
        with = FALSE
      ]
    )

    safe_write_csv(
      final_lock,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LOCKED_ARCHITECTURE_DEVELOPMENT_SPECIFICATION.csv"
        )
      )
    )

    safe_save_rds(
      final_genes,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_LOCKED_ARCHITECTURE_TOP50_GENES.rds"
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 18. STAGE-03 METHOD MANIFEST
# ==============================================================================

method_manifest <- data.table(
  experiment = c(
    "Experiment 1",
    "Experiment 2",
    "Experiment 3",
    "Experiment 4"
  ),
  development_design = c(
    "10-fold x 2-repeat CV",
    "5 outer x 10 inner nested CV; 20 MBO evaluations",
    "Development-only LOCO optimisation",
    "5 outer x 3 inner architecture CV"
  ),
  optimisation_target = c(
    "Fixed learner benchmarking",
    "Balanced accuracy",
    "Worst-cohort AUROC hierarchy",
    "Mean inner AUROC for LightGBM grid"
  ),
  held_out_final_evaluation = FALSE
)

safe_write_csv(
  method_manifest,
  file.path(
    AUDIT_ROOT,
    "TRACE_Omics_stage03_method_manifest.csv"
  )
)


# ==============================================================================
# 19. RUN REQUESTED EXPERIMENTS
# ==============================================================================

if (run_requested(
  "experiment1"
)) {
  run_experiment_1()
}

if (run_requested(
  "experiment2"
)) {
  run_experiment_2()
}

if (run_requested(
  "experiment3"
)) {
  run_experiment_3()
}

if (run_requested(
  "experiment4"
)) {
  run_experiment_4()
}


# ==============================================================================
# 20. FINAL REPRODUCIBILITY RECORDS
# ==============================================================================

config_snapshot <- file.path(
  SESSION_ROOT,
  paste0(
    "config_snapshot_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    ),
    ".yml"
  )
)

file.copy(
  CONFIG_PATH,
  config_snapshot,
  overwrite = TRUE
)

write_session_info(
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage03_model_development_optimisation_sessionInfo.txt"
  )
)

run_manifest <- data.table(
  timestamp = format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  config = CONFIG_PATH,
  modes = paste(
    RUN_MODES,
    collapse = ","
  ),
  diseases = paste(
    REQUESTED_DISEASES,
    collapse = ","
  ),
  seed = GLOBAL_SEED,
  experiment1_folds = 10L,
  experiment1_repeats = 2L,
  experiment2_outer_folds = CROSS_DISEASE_OUTER_FOLDS,
  experiment2_inner_folds = CROSS_DISEASE_INNER_FOLDS,
  experiment2_mbo_evaluations = CROSS_DISEASE_TUNING_EVALS,
  experiment3_feature_counts = paste(
    EXPERIMENT3_FEATURE_COUNTS,
    collapse = "/"
  ),
  experiment4_outer_folds = 5L,
  experiment4_inner_folds = 3L,
  experiment4_lightgbm_configurations = nrow(
    LIGHTGBM_GRID
  ),
  strict_audits = STRICT_AUDITS
)

safe_write_csv(
  run_manifest,
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage03_run_manifest.csv"
  )
)

cat(
  "\n============================================================\n",
  "TRACE-OMICS STAGE 03 COMPLETE\n",
  "Model development and optimisation finished successfully.\n",
  "Frozen final held-out evaluation remains reserved for Stage 04.\n",
  "Output root: ",
  STAGE3_ROOT,
  "\n============================================================\n",
  sep = ""
)

# ==============================================================================
# End of scripts/03_model_development_optimisation/
# 03_model_development_optimisation.R
# ==============================================================================
