# ==============================================================================
# TRACE-Omics
# Stage 04: Model testing and reliability assessment
# ==============================================================================
#
# File:
#   scripts/04_model_testing_reliability/04_model_testing_reliability.R
#
# Framework:
#   Transcriptomic Robustness Assessment and Cross-cohort Evaluation for
#   Omics learning (TRACE-Omics)
#
# Purpose:
#   Consume the development-only analytical locks produced by Stage 03 and
#   perform the testing, transportability and molecular-reliability analyses
#   described in the TRACE-Omics dissertation.
#
#   Experiment 1 - sample-support perturbation
#     * Fit each retained representation from the complete 66-sample
#       development population.
#     * Recalculate development-only RF permutation importance (300 trees).
#     * Freeze the top 50 genes and deterministic LightGBM classifier.
#     * Predict the untouched 16-sample balanced held-out population.
#     * Enumerate all choose(16, 8) = 12,870 balanced label allocations.
#     * Apply exact one-sided permutation testing and BH adjustment.
#     * Compare representations with paired DeLong tests.
#     * Quantify the development-to-held-out AUROC gap by 5,000 bootstrap
#       replicates using the actual Stage-03 development-fold results available.
#
#   Experiment 2 - disease-composition perturbation
#     * Reconstruct strict leave-one-disease-out transfer from source data only.
#     * Refit the locked Model A (TB-naive) and Model B (ALL5) systems from
#       their complete development populations.
#     * Evaluate both models on the identical frozen TB population.
#     * Compare Model A and Model B with paired DeLong testing.
#     * Estimate pointwise ROC uncertainty with 2,000 stratified bootstraps.
#     * Quantify gene-level Delta-AUROC importance with 50 permutations/gene.
#     * Contextualise positive predictive genes with Reactome enrichment using
#       the common 4,450-gene background.
#
#   Experiment 3 - cohort-composition perturbation
#     * Re-score the Stage-03 development-only LOCO predictions using the
#       frozen disease-specific operating threshold.
#     * Produce cohort-level ROC, PR and calibration evidence.
#     * Reconstruct each final locked Random Forest from the complete
#       development population and evaluate the untouched held-out population.
#     * Quantify final model reliance with 100 repeated Delta-AUROC
#       permutations per selected gene.
#     * Produce first-order ALE profiles for up to six leading genes.
#
#   Experiment 4 - model-architecture perturbation
#     * Evaluate LightGBM and TabPFN under identical molecular inputs.
#     * Reconstruct source-only architecture LOCO with source-only LightGBM
#       configuration selection.
#     * Fit the final locked LightGBM configuration and default pretrained
#       TabPFN to the same complete development matrix and the same 50 genes.
#     * Predict the same frozen diagnostic samples.
#     * Serialise and restore both Python model objects, requiring probability
#       reproduction within 1e-5 before interpretation.
#     * Quantify architecture-specific Delta-AUROC importance, Spearman rank
#       agreement, top-20 Jaccard overlap and model-specific ALE profiles.
#
# Critical stage boundary:
#   All representation, feature, class-handling, learner, hyperparameter and
#   operating-threshold choices are loaded from Stage 03 or reconstructed from
#   development-only evidence before the final held-out labels are scored.
#   Stage 04 never uses held-out performance to revise a predictive system.
#
# Recommended architecture execution:
#   Run Experiment 4 one disease at a time from a fresh R session so reticulate
#   can bind to the intended disease-specific Python environment.
#
# Example:
#
#   Rscript scripts/04_model_testing_reliability/04_model_testing_reliability.R \
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
      setdiff(RUN_MODES, VALID_MODES),
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
    "ranger",
    "pROC",
    "PRROC",
    "ggplot2"
  ),
  bioc = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

config <- read_trace_config(
  CONFIG_PATH
)

GLOBAL_SEED <- as.integer(
  get_config_value(
    config,
    "reproducibility.global_seed",
    get_config_value(
      config,
      "project.seed",
      7832
    )
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
      setdiff(REQUESTED_DISEASES, ALL_DISEASES),
      collapse = ", "
    ),
    call. = FALSE
  )
}

set_trace_seed(
  GLOBAL_SEED
)


# ==============================================================================
# 2. IMPORT STAGE-03 FUNCTION DEFINITIONS ONLY
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

STAGE3_SCRIPT <- file.path(
  "scripts",
  "03_model_development_optimisation",
  "03_model_development_optimisation.R"
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

  temporary_file <- tempfile(
    fileext = ".R"
  )

  on.exit(
    unlink(temporary_file),
    add = TRUE
  )

  writeLines(
    lines[seq_len(marker_line[[1]] - 1L)],
    temporary_file,
    useBytes = TRUE
  )

  sys.source(
    temporary_file,
    envir = environment
  )

  invisible(environment)
}

stage3 <- new.env(
  parent = globalenv()
)

source_prefix_before_marker(
  STAGE3_SCRIPT,
  "# 18. STAGE-03 METHOD MANIFEST",
  stage3
)

stage1 <- stage3$stage1
stage2 <- stage3$stage2

required_stage3_functions <- c(
  "load_canonical_disease",
  "load_equal_size_development",
  "classification_metrics",
  "safe_auroc",
  "safe_auprc",
  "make_fixed_learner",
  "make_mlr3_task",
  "prepare_cross_disease_development",
  "load_cross_disease_common_genes",
  "fit_predict_disease_rf",
  "class_weight_vector",
  "initialise_architecture_python"
)

for (function_name in required_stage3_functions) {
  if (!exists(
    function_name,
    envir = stage3,
    inherits = FALSE
  )) {
    stop(
      "Stage 03 function unavailable after import: ",
      function_name,
      call. = FALSE
    )
  }
}


# ==============================================================================
# 3. OUTPUT ROOTS AND REPRODUCIBILITY AUDITS
# ==============================================================================

stage4_example <- get_config_value(
  config,
  "paths.stage_results.generalisation_reliability",
  "results/04_generalisation_reliability/TB"
)

STAGE4_ROOT <- dirname(
  stage4_example
)

STAGE3_ROOT <- dirname(
  get_config_value(
    config,
    "paths.stage_results.model_development_optimisation",
    "results/03_model_development_optimisation/TB"
  )
)

AUDIT_ROOT <- file.path(
  STAGE4_ROOT,
  "00_audits"
)

SESSION_ROOT <- file.path(
  STAGE4_ROOT,
  "99_session_info"
)

for (path in c(
  STAGE4_ROOT,
  AUDIT_ROOT,
  SESSION_ROOT
)) {
  make_dir(path)
}

prerequisite_paths <- c(
  STAGE1_SCRIPT,
  STAGE2_SCRIPT,
  STAGE3_SCRIPT
)

stage_source_manifest <- data.table(
  stage = c(
    "Stage01",
    "Stage02",
    "Stage03"
  ),
  path = prerequisite_paths,
  md5 = unname(
    tools::md5sum(
      prerequisite_paths
    )
  )
)

safe_write_csv(
  stage_source_manifest,
  file.path(
    AUDIT_ROOT,
    "TRACE_Omics_stage04_prerequisite_md5.csv"
  )
)

message_section(
  "TRACE-Omics Stage 04: Model testing and reliability assessment"
)

message(
  "Configuration: ",
  CONFIG_PATH
)

message(
  "Modes:         ",
  paste(RUN_MODES, collapse = ", ")
)

message(
  "Diseases:      ",
  paste(REQUESTED_DISEASES, collapse = ", ")
)

message(
  "Global seed:   ",
  GLOBAL_SEED
)


# ==============================================================================
# 4. FROZEN-PARTITION LOADERS AND AUDITS
# ==============================================================================

normalise_binary_metadata <- function(metadata) {
  metadata <- as.data.frame(
    metadata,
    check.names = FALSE
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

  metadata
}

subset_partition_from_manifest <- function(canonical,
                                           role) {
  ids <- as.character(
    canonical$manifest$sample_id[
      canonical$manifest$role == role
    ]
  )

  metadata <- canonical$metadata[
    match(
      ids,
      canonical$metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata$sample_id)) {
    stop(
      canonical$disease,
      " ",
      role,
      " manifest contains IDs absent from canonical metadata.",
      call. = FALSE
    )
  }

  metadata <- normalise_binary_metadata(
    metadata
  )

  counts <- canonical$counts[
    ,
    metadata$sample_id,
    drop = FALSE
  ]

  if (!identical(
    colnames(counts),
    metadata$sample_id
  )) {
    stop(
      canonical$disease,
      " ",
      role,
      " sample alignment failed.",
      call. = FALSE
    )
  }

  list(
    counts = counts,
    metadata = metadata
  )
}

load_general_split <- function(disease_id) {
  canonical <- stage3$load_canonical_disease(
    disease_id
  )

  development <- subset_partition_from_manifest(
    canonical,
    "development"
  )

  held_out <- subset_partition_from_manifest(
    canonical,
    "held_out"
  )

  overlap <- intersect(
    development$metadata$sample_id,
    held_out$metadata$sample_id
  )

  audit <- data.table(
    disease = disease_id,
    development_n = nrow(development$metadata),
    held_out_n = nrow(held_out$metadata),
    total_manifest_n = nrow(canonical$manifest),
    total_canonical_n = nrow(canonical$metadata),
    overlap_n = length(overlap),
    complete_accounting =
      nrow(development$metadata) +
        nrow(held_out$metadata) ==
        nrow(canonical$metadata)
  )

  if (
    audit$overlap_n[[1]] != 0L ||
      !audit$complete_accounting[[1]]
  ) {
    stop(
      disease_id,
      " general development/held-out partition audit failed.",
      call. = FALSE
    )
  }

  list(
    canonical = canonical,
    development = development,
    held_out = held_out,
    audit = audit
  )
}

load_equal_size_split <- function(disease_id) {
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

  canonical <- stage3$load_canonical_disease(
    disease_id
  )

  manifest <- as.data.frame(
    safe_read_csv(
      manifest_path
    ),
    check.names = FALSE
  )

  check_required_columns(
    manifest,
    c(
      "sample_id",
      "role"
    ),
    paste0(
      disease_id,
      " equal-size manifest"
    )
  )

  make_part <- function(role) {
    ids <- as.character(
      manifest$sample_id[
        manifest$role == role
      ]
    )

    metadata <- canonical$metadata[
      match(
        ids,
        canonical$metadata$sample_id
      ),
      ,
      drop = FALSE
    ]

    metadata <- normalise_binary_metadata(
      metadata
    )

    counts <- canonical$counts[
      ,
      metadata$sample_id,
      drop = FALSE
    ]

    list(
      counts = counts,
      metadata = metadata
    )
  }

  development <- make_part(
    "development"
  )

  held_out <- make_part(
    "held_out"
  )

  audit <- data.table(
    disease = disease_id,
    development_n = nrow(development$metadata),
    development_cases = sum(
      development$metadata$outcome ==
        POSITIVE_CLASS
    ),
    development_controls = sum(
      development$metadata$outcome ==
        NEGATIVE_CLASS
    ),
    held_out_n = nrow(held_out$metadata),
    held_out_cases = sum(
      held_out$metadata$outcome ==
        POSITIVE_CLASS
    ),
    held_out_controls = sum(
      held_out$metadata$outcome ==
        NEGATIVE_CLASS
    ),
    overlap_n = length(
      intersect(
        development$metadata$sample_id,
        held_out$metadata$sample_id
      )
    )
  )

  expected_ok <- (
    audit$development_n[[1]] == 66L &&
      audit$development_cases[[1]] == 33L &&
      audit$development_controls[[1]] == 33L &&
      audit$held_out_n[[1]] == 16L &&
      audit$held_out_cases[[1]] == 8L &&
      audit$held_out_controls[[1]] == 8L &&
      audit$overlap_n[[1]] == 0L
  )

  audit[
    ,
    expected_design_ok := expected_ok
  ]

  if (!expected_ok) {
    stop(
      disease_id,
      " equal-size split differs from the prespecified 66/16 balanced design.",
      call. = FALSE
    )
  }

  list(
    canonical = canonical,
    manifest = manifest,
    development = development,
    held_out = held_out,
    audit = audit
  )
}

audit_expected_cross_disease_counts <- function(model_a_n,
                                                model_b_n,
                                                tb_test_n) {
  expected_a <- as.integer(
    get_config_value(
      config,
      "disease_composition_experiment.model_a.expected_development_n",
      1413L
    )
  )

  expected_b <- as.integer(
    get_config_value(
      config,
      "disease_composition_experiment.model_b.expected_development_n",
      1736L
    )
  )

  audit <- data.table(
    quantity = c(
      "Model_A_development_n",
      "Model_B_development_n",
      "Frozen_TB_test_n"
    ),
    observed = c(
      model_a_n,
      model_b_n,
      tb_test_n
    ),
    expected = c(
      expected_a,
      expected_b,
      142L
    )
  )

  audit[
    ,
    match := observed == expected
  ]

  safe_write_csv(
    audit,
    file.path(
      AUDIT_ROOT,
      "Experiment2_cross_disease_sample_count_audit.csv"
    )
  )

  if (!all(audit$match)) {
    warning(
      "Experiment 2 sample counts differ from the dissertation reference values. ",
      "Stage 04 will retain the frozen manifests generated by Stage 01 rather ",
      "than silently construct a second split."
    )

    if (STRICT_AUDITS) {
      stop(
        "Strict audit requested and Experiment 2 sample-count reference did not match.",
        call. = FALSE
      )
    }
  }

  audit
}


# ==============================================================================
# 5. PERFORMANCE, CURVE AND CALIBRATION HELPERS
# ==============================================================================

truth_numeric <- function(truth) {
  as.integer(
    as.character(truth) ==
      POSITIVE_CLASS
  )
}

safe_delong_ci <- function(truth,
                           probability) {
  truth <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  if (
    length(unique(truth)) < 2L ||
      length(probability) != length(truth)
  ) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

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

  interval <- as.numeric(
    pROC::ci.auc(
      roc_object,
      method = "delong"
    )
  )

  c(
    lower = interval[[1]],
    upper = interval[[3]]
  )
}

performance_metrics <- function(truth,
                                probability,
                                threshold = 0.50) {
  base <- stage3$classification_metrics(
    truth = truth,
    probability = probability,
    threshold = threshold
  )

  ci <- safe_delong_ci(
    truth,
    probability
  )

  base[
    ,
    `:=`(
      AUROC_lower95 = ci[["lower"]],
      AUROC_upper95 = ci[["upper"]],
      Prevalence = mean(
        as.character(truth) ==
          POSITIVE_CLASS
      ),
      Threshold = threshold
    )
  ]

  base
}

empirical_roc_curve <- function(truth,
                                probability) {
  roc_object <- pROC::roc(
    response = factor(
      as.character(truth),
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    predictor = probability,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  data.table(
    specificity = as.numeric(
      roc_object$specificities
    ),
    sensitivity = as.numeric(
      roc_object$sensitivities
    )
  )[
    ,
    `:=`(
      FPR = 1 - specificity,
      TPR = sensitivity
    )
  ][
    order(FPR, TPR)
  ]
}

empirical_pr_curve <- function(truth,
                               probability) {
  y <- truth_numeric(
    truth
  )

  ord <- order(
    probability,
    decreasing = TRUE
  )

  y <- y[ord]

  tp <- cumsum(
    y == 1L
  )

  fp <- cumsum(
    y == 0L
  )

  positives <- sum(
    y == 1L
  )

  recall <- if (positives > 0L) {
    tp / positives
  } else {
    rep(
      NA_real_,
      length(tp)
    )
  }

  precision <- tp /
    pmax(
      tp + fp,
      1L
    )

  data.table(
    recall = c(
      0,
      recall
    ),
    precision = c(
      1,
      precision
    )
  )
}

make_calibration_data <- function(truth,
                                  probability,
                                  bins = 6L) {
  y <- truth_numeric(
    truth
  )

  dt <- data.table(
    truth = y,
    probability = as.numeric(
      probability
    )
  )

  dt[
    ,
    rank_index := data.table::frank(
      probability,
      ties.method = "average"
    )
  ]

  dt[
    ,
    calibration_bin := pmin(
      bins,
      pmax(
        1L,
        ceiling(
          rank_index /
            .N *
            bins
        )
      )
    )
  ]

  output <- dt[
    ,
    {
      successes <- sum(
        truth == 1L
      )

      n_bin <- .N

      interval <- stats::binom.test(
        successes,
        n_bin
      )$conf.int

      .(
        n = n_bin,
        mean_predicted = mean(probability),
        observed_fraction = mean(truth),
        lower95 = interval[[1]],
        upper95 = interval[[2]]
      )
    },
    by = calibration_bin
  ]

  setorder(
    output,
    calibration_bin
  )

  output
}

calibration_intercept_slope <- function(truth,
                                        probability) {
  y <- truth_numeric(
    truth
  )

  epsilon <- 1e-6

  p <- pmin(
    1 - epsilon,
    pmax(
      epsilon,
      as.numeric(probability)
    )
  )

  lp <- stats::qlogis(
    p
  )

  slope <- tryCatch(
    {
      model <- stats::glm(
        y ~ lp,
        family = stats::binomial()
      )

      unname(
        stats::coef(model)[["lp"]]
      )
    },
    error = function(e) NA_real_
  )

  intercept <- tryCatch(
    {
      model <- stats::glm(
        y ~ 1,
        offset = lp,
        family = stats::binomial()
      )

      unname(
        stats::coef(model)[["(Intercept)"]]
      )
    },
    error = function(e) NA_real_
  )

  data.table(
    Calibration_intercept = intercept,
    Calibration_slope = slope
  )
}

stratified_resample_indices <- function(truth) {
  truth_chr <- as.character(
    truth
  )

  case_index <- which(
    truth_chr ==
      POSITIVE_CLASS
  )

  control_index <- which(
    truth_chr ==
      NEGATIVE_CLASS
  )

  c(
    sample(
      case_index,
      length(case_index),
      replace = TRUE
    ),
    sample(
      control_index,
      length(control_index),
      replace = TRUE
    )
  )
}

interpolate_step <- function(x,
                             y,
                             grid) {
  if (
    length(x) < 2L ||
      all(!is.finite(x)) ||
      all(!is.finite(y))
  ) {
    return(
      rep(
        NA_real_,
        length(grid)
      )
    )
  }

  x <- as.numeric(x)
  y <- as.numeric(y)

  keep <- is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  ord <- order(
    x,
    y
  )

  x <- x[ord]
  y <- y[ord]

  unique_x <- unique(
    x
  )

  aggregate_y <- vapply(
    unique_x,
    function(current_x) {
      max(
        y[x == current_x],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  stats::approx(
    x = unique_x,
    y = aggregate_y,
    xout = grid,
    rule = 2,
    ties = "ordered"
  )$y
}

bootstrap_curve_bands <- function(truth,
                                  probability,
                                  iterations = 2000L,
                                  seed = GLOBAL_SEED,
                                  grid_points = 101L) {
  set.seed(
    seed
  )

  specificity_grid <- seq(
    0,
    1,
    length.out = grid_points
  )

  recall_grid <- seq(
    0,
    1,
    length.out = grid_points
  )

  roc_matrix <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = grid_points
  )

  pr_matrix <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = grid_points
  )

  for (iteration in seq_len(iterations)) {
    index <- stratified_resample_indices(
      truth
    )

    boot_truth <- truth[index]
    boot_probability <- probability[index]

    roc_dt <- empirical_roc_curve(
      boot_truth,
      boot_probability
    )

    roc_matrix[
      iteration,
    ] <- interpolate_step(
      x = roc_dt$specificity,
      y = roc_dt$sensitivity,
      grid = specificity_grid
    )

    pr_dt <- empirical_pr_curve(
      boot_truth,
      boot_probability
    )

    pr_matrix[
      iteration,
    ] <- interpolate_step(
      x = pr_dt$recall,
      y = pr_dt$precision,
      grid = recall_grid
    )
  }

  quantile_column <- function(x,
                              probability_level) {
    apply(
      x,
      2L,
      stats::quantile,
      probs = probability_level,
      na.rm = TRUE,
      names = FALSE
    )
  }

  list(
    roc = data.table(
      specificity = specificity_grid,
      lower95 = quantile_column(
        roc_matrix,
        0.025
      ),
      upper95 = quantile_column(
        roc_matrix,
        0.975
      )
    )[
      ,
      FPR := 1 - specificity
    ],
    pr = data.table(
      recall = recall_grid,
      lower95 = quantile_column(
        pr_matrix,
        0.025
      ),
      upper95 = quantile_column(
        pr_matrix,
        0.975
      )
    )
  )
}

paired_delong_test <- function(truth,
                               probability_a,
                               probability_b,
                               label_a,
                               label_b) {
  truth_factor <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  roc_a <- pROC::roc(
    truth_factor,
    probability_a,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  roc_b <- pROC::roc(
    truth_factor,
    probability_b,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  test <- pROC::roc.test(
    roc_a,
    roc_b,
    paired = TRUE,
    method = "delong",
    alternative = "two.sided"
  )

  data.table(
    Model_A = label_a,
    Model_B = label_b,
    AUROC_A = as.numeric(
      pROC::auc(roc_a)
    ),
    AUROC_B = as.numeric(
      pROC::auc(roc_b)
    ),
    Delta_A_minus_B =
      as.numeric(pROC::auc(roc_a)) -
        as.numeric(pROC::auc(roc_b)),
    p_value = as.numeric(
      test$p.value
    )
  )
}

save_plot_dual <- function(plot,
                           path_stem,
                           width = 7,
                           height = 6) {
  make_dir(
    dirname(path_stem)
  )

  ggplot2::ggsave(
    filename = paste0(
      path_stem,
      ".pdf"
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggplot2::ggsave(
    filename = paste0(
      path_stem,
      ".png"
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}

save_performance_bundle <- function(output_dir,
                                    prefix,
                                    truth,
                                    probability,
                                    threshold = 0.50,
                                    calibration_bins = 6L,
                                    bootstrap_iterations = 2000L,
                                    seed = GLOBAL_SEED,
                                    save_bootstrap = TRUE) {
  make_dir(
    output_dir
  )

  metrics <- performance_metrics(
    truth,
    probability,
    threshold
  )

  calibration_stats <- calibration_intercept_slope(
    truth,
    probability
  )

  metrics <- cbind(
    metrics,
    calibration_stats
  )

  roc_dt <- empirical_roc_curve(
    truth,
    probability
  )

  pr_dt <- empirical_pr_curve(
    truth,
    probability
  )

  calibration_dt <- make_calibration_data(
    truth,
    probability,
    bins = calibration_bins
  )

  safe_write_csv(
    metrics,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_metrics.csv"
      )
    )
  )

  safe_write_csv(
    roc_dt,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_ROC_coordinates.csv"
      )
    )
  )

  safe_write_csv(
    pr_dt,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_PR_coordinates.csv"
      )
    )
  )

  safe_write_csv(
    calibration_dt,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_calibration_coordinates.csv"
      )
    )
  )

  bands <- NULL

  if (
    isTRUE(save_bootstrap) &&
      bootstrap_iterations > 0L
  ) {
    bands <- bootstrap_curve_bands(
      truth = truth,
      probability = probability,
      iterations = bootstrap_iterations,
      seed = seed,
      grid_points = as.integer(
        get_config_value(
          config,
          "performance.curve_uncertainty.grid_points",
          101L
        )
      )
    )

    safe_write_csv(
      bands$roc,
      file.path(
        output_dir,
        paste0(
          prefix,
          "_ROC_bootstrap_band.csv"
        )
      )
    )

    safe_write_csv(
      bands$pr,
      file.path(
        output_dir,
        paste0(
          prefix,
          "_PR_bootstrap_band.csv"
        )
      )
    )
  }

  prevalence <- mean(
    as.character(truth) ==
      POSITIVE_CLASS
  )

  roc_plot <- ggplot2::ggplot(
    roc_dt,
    ggplot2::aes(
      x = FPR,
      y = TPR
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_step(
      linewidth = 0.8
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      title = prefix
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  pr_plot <- ggplot2::ggplot(
    pr_dt,
    ggplot2::aes(
      x = recall,
      y = precision
    )
  ) +
    ggplot2::geom_hline(
      yintercept = prevalence,
      linetype = 2
    ) +
    ggplot2::geom_step(
      linewidth = 0.8
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Recall",
      y = "Precision",
      title = prefix
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  calibration_plot <- ggplot2::ggplot(
    calibration_dt,
    ggplot2::aes(
      x = mean_predicted,
      y = observed_fraction
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower95,
        ymax = upper95
      ),
      width = 0
    ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Mean predicted probability",
      y = "Observed case fraction",
      title = prefix
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  figure_dir <- file.path(
    output_dir,
    "figures"
  )

  save_plot_dual(
    roc_plot,
    file.path(
      figure_dir,
      paste0(
        prefix,
        "_ROC"
      )
    )
  )

  save_plot_dual(
    pr_plot,
    file.path(
      figure_dir,
      paste0(
        prefix,
        "_PR"
      )
    )
  )

  save_plot_dual(
    calibration_plot,
    file.path(
      figure_dir,
      paste0(
        prefix,
        "_calibration"
      )
    )
  )

  list(
    metrics = metrics,
    roc = roc_dt,
    pr = pr_dt,
    calibration = calibration_dt,
    bands = bands
  )
}



# ==============================================================================
# 6. FITTED-MODEL AND MOLECULAR-RELIABILITY HELPERS
# ==============================================================================

fit_fixed_mlr3_model <- function(learner_id,
                                 train_matrix,
                                 train_metadata,
                                 seed) {
  task_info <- stage3$make_mlr3_task(
    expression_matrix = train_matrix,
    metadata = train_metadata,
    task_id = paste0(
      "stage04_",
      safe_file_name(learner_id)
    )
  )

  learner <- stage3$make_fixed_learner(
    learner_id,
    seed
  )

  set.seed(
    seed
  )

  learner$train(
    task_info$task
  )

  list(
    learner = learner,
    gene_map = task_info$map
  )
}

predict_fixed_mlr3_model <- function(model_object,
                                     evaluation_matrix) {
  gene_map <- model_object$gene_map

  missing_genes <- setdiff(
    gene_map$original_gene,
    rownames(evaluation_matrix)
  )

  if (length(missing_genes) > 0L) {
    stop(
      "Evaluation matrix is missing fitted genes: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  x <- as.data.frame(
    t(
      evaluation_matrix[
        gene_map$original_gene,
        ,
        drop = FALSE
      ]
    ),
    check.names = FALSE
  )

  colnames(x) <- gene_map$safe_gene

  prediction <- model_object$learner$predict_newdata(
    x
  )

  as.numeric(
    prediction$prob[
      ,
      POSITIVE_CLASS
    ]
  )
}

fit_rf_importance_ranking <- function(train_matrix,
                                      train_metadata,
                                      trees,
                                      top_n,
                                      seed) {
  stage2$select_ranked_features(
    train_matrix = train_matrix,
    train_metadata = train_metadata,
    filter_id = "importance",
    top_n = top_n,
    seed = seed,
    rf_trees = trees
  )
}

fit_final_ranger_model <- function(train_matrix,
                                   train_metadata,
                                   weighting,
                                   seed,
                                   trees = 250L,
                                   min_node_size = 5L) {
  train_matrix <- coerce_numeric_matrix(
    train_matrix,
    "final Random Forest training matrix"
  )

  x <- as.data.frame(
    t(train_matrix),
    check.names = FALSE
  )

  original_genes <- colnames(x)

  safe_genes <- make.unique(
    make.names(
      original_genes,
      unique = TRUE
    )
  )

  colnames(x) <- safe_genes

  x$outcome <- factor(
    as.character(
      train_metadata$outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  weights <- if (weighting == "Weighted") {
    stage3$class_weight_vector(
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
    data = x,
    probability = TRUE,
    num.trees = as.integer(trees),
    mtry = max(
      1L,
      floor(
        sqrt(
          nrow(train_matrix)
        )
      )
    ),
    min.node.size = as.integer(
      min_node_size
    ),
    class.weights = weights,
    importance = "none",
    num.threads = 1L,
    seed = as.integer(seed)
  )

  list(
    model = model,
    genes = rownames(train_matrix),
    safe_genes = safe_genes,
    weighting = weighting
  )
}

predict_final_ranger_model <- function(model_object,
                                       matrix) {
  missing <- setdiff(
    model_object$genes,
    rownames(matrix)
  )

  if (length(missing) > 0L) {
    stop(
      "Random Forest prediction matrix is missing fitted genes: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x <- as.data.frame(
    t(
      matrix[
        model_object$genes,
        ,
        drop = FALSE
      ]
    ),
    check.names = FALSE
  )

  colnames(x) <- model_object$safe_genes

  prediction <- predict(
    model_object$model,
    data = x
  )$predictions

  as.numeric(
    prediction[
      ,
      POSITIVE_CLASS
    ]
  )
}

fit_python_lightgbm_object <- function(x_train,
                                       y_train,
                                       parameters,
                                       seed,
                                       python) {
  model <- python$lightgbm$LGBMClassifier(
    objective = "binary",
    n_estimators = as.integer(
      parameters$n_estimators[[1L]]
    ),
    learning_rate = as.numeric(
      parameters$learning_rate[[1L]]
    ),
    num_leaves = as.integer(
      parameters$num_leaves[[1L]]
    ),
    max_depth = as.integer(
      parameters$max_depth[[1L]]
    ),
    min_child_samples = as.integer(
      parameters$min_child_samples[[1L]]
    ),
    subsample = as.numeric(
      parameters$subsample[[1L]]
    ),
    subsample_freq = 1L,
    colsample_bytree = as.numeric(
      parameters$colsample_bytree[[1L]]
    ),
    reg_alpha = as.numeric(
      parameters$reg_alpha[[1L]]
    ),
    reg_lambda = as.numeric(
      parameters$reg_lambda[[1L]]
    ),
    class_weight = "balanced",
    random_state = as.integer(seed),
    n_jobs = -1L,
    verbosity = -1L
  )

  model$fit(
    stage3$to_numpy_x(
      x_train,
      python
    ),
    stage3$to_numpy_y(
      y_train,
      python
    )
  )

  model
}

fit_python_tabpfn_object <- function(x_train,
                                     y_train,
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
    stage3$to_numpy_x(
      x_train,
      python
    ),
    stage3$to_numpy_y(
      y_train,
      python
    )
  )

  model
}

predict_python_object <- function(model,
                                  x_test,
                                  python) {
  probability_matrix <- reticulate::py_to_r(
    model$predict_proba(
      stage3$to_numpy_x(
        x_test,
        python
      )
    )
  )

  if (!is.matrix(probability_matrix)) {
    probability_matrix <- as.matrix(
      probability_matrix
    )
  }

  if (ncol(probability_matrix) != 2L) {
    stop(
      "Expected two probability columns from Python classifier; observed ",
      ncol(probability_matrix),
      ".",
      call. = FALSE
    )
  }

  as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )
}

make_python_prediction_function <- function(model,
                                            python) {
  force(model)
  force(python)

  function(model,
           newdata) {
    predict_python_object(
      model = model,
      x_test = as.matrix(newdata),
      python = python
    )
  }
}

make_ranger_prediction_function <- function(model_object) {
  force(model_object)

  function(model,
           newdata) {
    x <- as.data.frame(
      newdata,
      check.names = FALSE
    )

    colnames(x) <- model_object$safe_genes

    prediction <- predict(
      model_object$model,
      data = x
    )$predictions

    as.numeric(
      prediction[
        ,
        POSITIVE_CLASS
      ]
    )
  }
}

run_delta_auroc_permutation <- function(model_label,
                                        x,
                                        truth,
                                        predict_probability,
                                        permutations,
                                        seed,
                                        output_dir,
                                        resume = TRUE) {
  make_dir(
    output_dir
  )

  x <- as.data.frame(
    x,
    check.names = FALSE
  )

  genes <- colnames(x)

  if (
    is.null(genes) ||
      anyNA(genes) ||
      any(!nzchar(genes))
  ) {
    stop(
      model_label,
      ": predictor matrix requires non-empty gene names.",
      call. = FALSE
    )
  }

  baseline_probability <- predict_probability(
    x
  )

  baseline_auroc <- stage3$safe_auroc(
    truth,
    baseline_probability
  )

  raw_paths <- file.path(
    output_dir,
    paste0(
      safe_file_name(model_label),
      "__",
      safe_file_name(genes),
      "__raw.csv"
    )
  )

  completed_rows <- vector(
    "list",
    length(genes)
  )

  for (gene_index in seq_along(genes)) {
    gene <- genes[[gene_index]]
    raw_path <- raw_paths[[gene_index]]

    if (
      isTRUE(resume) &&
        file.exists(raw_path)
    ) {
      existing <- data.table::fread(
        raw_path
      )

      if (
        nrow(existing) == permutations &&
          all(existing$Gene == gene)
      ) {
        completed_rows[[gene_index]] <- existing
        next
      }
    }

    rows <- vector(
      "list",
      permutations
    )

    for (permutation_index in seq_len(permutations)) {
      current_seed <- make_trace_subseed(
        seed,
        model_label,
        gene,
        permutation_index
      )

      set.seed(
        current_seed
      )

      x_permuted <- x
      x_permuted[[gene]] <- sample(
        x_permuted[[gene]],
        replace = FALSE
      )

      permuted_probability <- predict_probability(
        x_permuted
      )

      permuted_auroc <- stage3$safe_auroc(
        truth,
        permuted_probability
      )

      rows[[permutation_index]] <- data.table(
        Model = model_label,
        Gene = gene,
        Permutation = permutation_index,
        Baseline_AUROC = baseline_auroc,
        Permuted_AUROC = permuted_auroc,
        Delta_AUROC =
          baseline_auroc -
            permuted_auroc,
        Seed = current_seed
      )
    }

    raw_dt <- data.table::rbindlist(
      rows,
      use.names = TRUE,
      fill = TRUE
    )

    safe_write_csv(
      raw_dt,
      raw_path
    )

    completed_rows[[gene_index]] <- raw_dt

    message(
      model_label,
      " permutation importance: ",
      gene_index,
      "/",
      length(genes),
      " genes complete."
    )
  }

  raw_dt <- data.table::rbindlist(
    completed_rows,
    use.names = TRUE,
    fill = TRUE
  )

  summary_dt <- raw_dt[
    ,
    .(
      Baseline_AUROC = unique(Baseline_AUROC)[[1L]],
      Mean_permuted_AUROC = mean(
        Permuted_AUROC,
        na.rm = TRUE
      ),
      Mean_Delta_AUROC = mean(
        Delta_AUROC,
        na.rm = TRUE
      ),
      SD_Delta_AUROC = stats::sd(
        Delta_AUROC,
        na.rm = TRUE
      ),
      Lower95_Delta_AUROC = as.numeric(
        stats::quantile(
          Delta_AUROC,
          probs = 0.025,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      Upper95_Delta_AUROC = as.numeric(
        stats::quantile(
          Delta_AUROC,
          probs = 0.975,
          na.rm = TRUE,
          names = FALSE
        )
      ),
      Positive_permutation_fraction = mean(
        Delta_AUROC > 0,
        na.rm = TRUE
      )
    ),
    by = .(
      Model,
      Gene
    )
  ]

  setorder(
    summary_dt,
    -Mean_Delta_AUROC,
    Gene
  )

  summary_dt[
    ,
    Rank := seq_len(.N)
  ]

  safe_write_csv(
    raw_dt,
    file.path(
      output_dir,
      paste0(
        safe_file_name(model_label),
        "_permutation_importance_raw.csv"
      )
    )
  )

  safe_write_csv(
    summary_dt,
    file.path(
      output_dir,
      paste0(
        safe_file_name(model_label),
        "_permutation_importance_summary.csv"
      )
    )
  )

  list(
    baseline_probability = baseline_probability,
    baseline_auroc = baseline_auroc,
    raw = raw_dt,
    summary = summary_dt
  )
}

make_dalex_explainer <- function(model,
                                 x,
                                 truth,
                                 predict_function,
                                 label) {
  require_packages(
    c(
      "DALEX",
      "ingredients"
    ),
    bioc = FALSE
  )

  DALEX::explain(
    model = model,
    data = as.data.frame(
      x,
      check.names = FALSE
    ),
    y = truth_numeric(
      truth
    ),
    predict_function = predict_function,
    label = label,
    verbose = FALSE
  )
}

run_ale_profiles <- function(explainer,
                             genes,
                             output_dir,
                             prefix) {
  require_packages(
    c(
      "DALEX",
      "ingredients"
    ),
    bioc = FALSE
  )

  make_dir(
    output_dir
  )

  genes <- unique(
    as.character(genes)
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  if (length(genes) == 0L) {
    warning(
      prefix,
      ": no genes were supplied for ALE."
    )

    return(
      invisible(NULL)
    )
  }

  set.seed(
    GLOBAL_SEED
  )

  ale <- ingredients::accumulated_dependency(
    explainer,
    variables = genes,
    N = NULL
  )

  safe_save_rds(
    ale,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_ALE.rds"
      )
    )
  )

  ale_dt <- data.table::as.data.table(
    ale
  )

  safe_write_csv(
    ale_dt,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_ALE.csv"
      )
    )
  )

  if (
    all(
      c(
        "_x_",
        "_yhat_",
        "_vname_"
      ) %in% names(ale_dt)
    )
  ) {
    plot <- ggplot2::ggplot(
      ale_dt,
      ggplot2::aes(
        x = `_x_`,
        y = `_yhat_`
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = 2
      ) +
      ggplot2::geom_line(
        linewidth = 0.8
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(`_vname_`),
        scales = "free_x"
      ) +
      ggplot2::labs(
        x = "Observed gene expression",
        y = "Accumulated local effect",
        title = prefix
      ) +
      ggplot2::theme_classic(
        base_size = 11
      )

    save_plot_dual(
      plot,
      file.path(
        output_dir,
        paste0(
          prefix,
          "_ALE"
        )
      ),
      width = 9,
      height = 7
    )
  }

  ale
}

run_reactome_enrichment <- function(importance_summary,
                                    background_genes,
                                    output_dir,
                                    prefix,
                                    top_positive_genes = 25L) {
  make_dir(
    output_dir
  )

  required <- c(
    "AnnotationDbi",
    "org.Hs.eg.db",
    "ReactomePA"
  )

  available <- vapply(
    required,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )

  package_audit <- data.table(
    package = required,
    available = available
  )

  safe_write_csv(
    package_audit,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_Reactome_package_audit.csv"
      )
    )
  )

  if (!all(available)) {
    warning(
      prefix,
      ": Reactome enrichment skipped because required Bioconductor packages are unavailable."
    )

    safe_write_csv(
      data.table(
        status = "skipped",
        reason = "Required Bioconductor package unavailable"
      ),
      file.path(
        output_dir,
        paste0(
          prefix,
          "_Reactome_status.csv"
        )
      )
    )

    return(
      invisible(NULL)
    )
  }

  selected_symbols <- importance_summary[
    is.finite(Mean_Delta_AUROC) &
      Mean_Delta_AUROC > 0
  ][
    order(-Mean_Delta_AUROC)
  ][
    seq_len(
      min(
        .N,
        top_positive_genes
      )
    ),
    Gene
  ]

  if (length(selected_symbols) == 0L) {
    safe_write_csv(
      data.table(
        status = "skipped",
        reason = "No genes had positive mean Delta-AUROC"
      ),
      file.path(
        output_dir,
        paste0(
          prefix,
          "_Reactome_status.csv"
        )
      )
    )

    return(
      invisible(NULL)
    )
  }

  selected_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(selected_symbols),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )

  background_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(background_genes),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )

  selected_entrez <- unique(
    selected_map$ENTREZID[
      !is.na(selected_map$ENTREZID)
    ]
  )

  background_entrez <- unique(
    background_map$ENTREZID[
      !is.na(background_map$ENTREZID)
    ]
  )

  if (
    length(selected_entrez) == 0L ||
      length(background_entrez) == 0L
  ) {
    safe_write_csv(
      data.table(
        status = "skipped",
        reason = "Gene-symbol to Entrez mapping yielded no usable identifiers"
      ),
      file.path(
        output_dir,
        paste0(
          prefix,
          "_Reactome_status.csv"
        )
      )
    )

    return(
      invisible(NULL)
    )
  }

  enrichment <- ReactomePA::enrichPathway(
    gene = selected_entrez,
    organism = "human",
    universe = background_entrez,
    minGSSize = 5,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    readable = TRUE
  )

  enrichment_dt <- data.table::as.data.table(
    as.data.frame(
      enrichment
    )
  )

  safe_write_csv(
    enrichment_dt,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_Reactome_enrichment.csv"
      )
    )
  )

  enrichment_dt
}


# ==============================================================================
# 7. EXPERIMENT 1 - SAMPLE-SUPPORT FINAL TESTING
# ==============================================================================

exact_rank_auroc <- function(truth_binary,
                             probability) {
  truth_binary <- as.integer(
    truth_binary
  )

  n_positive <- sum(
    truth_binary == 1L
  )

  n_negative <- sum(
    truth_binary == 0L
  )

  if (
    n_positive == 0L ||
      n_negative == 0L
  ) {
    return(
      NA_real_
    )
  }

  ranks <- rank(
    probability,
    ties.method = "average"
  )

  (
    sum(
      ranks[truth_binary == 1L]
    ) -
      n_positive *
        (n_positive + 1) /
        2
  ) /
    (
      n_positive *
        n_negative
    )
}

exact_balanced_permutation_test <- function(probability,
                                            observed_truth) {
  y <- truth_numeric(
    observed_truth
  )

  n <- length(
    y
  )

  n_positive <- sum(
    y == 1L
  )

  combinations <- utils::combn(
    n,
    n_positive
  )

  observed_auroc <- exact_rank_auroc(
    y,
    probability
  )

  null_auroc <- vapply(
    seq_len(
      ncol(combinations)
    ),
    function(index) {
      permuted <- integer(
        n
      )

      permuted[
        combinations[
          ,
          index
        ]
      ] <- 1L

      exact_rank_auroc(
        permuted,
        probability
      )
    },
    numeric(1)
  )

  data.table(
    observed_AUROC = observed_auroc,
    null_allocation_count = length(null_auroc),
    exact_one_sided_p = mean(
      null_auroc >=
        observed_auroc
    )
  ) -> result

  list(
    result = result,
    null = data.table(
      allocation = seq_along(null_auroc),
      null_AUROC = null_auroc
    )
  )
}

bootstrap_development_to_heldout_gap <- function(development_fold_aurocs,
                                                  heldout_truth,
                                                  heldout_probability,
                                                  iterations = 5000L,
                                                  seed = GLOBAL_SEED) {
  development_fold_aurocs <- as.numeric(
    development_fold_aurocs
  )

  development_fold_aurocs <- development_fold_aurocs[
    is.finite(
      development_fold_aurocs
    )
  ]

  if (length(development_fold_aurocs) == 0L) {
    stop(
      "Development-to-held-out bootstrap requires at least one finite development-fold AUROC.",
      call. = FALSE
    )
  }

  set.seed(
    seed
  )

  gap <- numeric(
    iterations
  )

  train_auc <- numeric(
    iterations
  )

  heldout_auc <- numeric(
    iterations
  )

  for (iteration in seq_len(iterations)) {
    train_auc[[iteration]] <- mean(
      sample(
        development_fold_aurocs,
        length(development_fold_aurocs),
        replace = TRUE
      )
    )

    index <- stratified_resample_indices(
      heldout_truth
    )

    heldout_auc[[iteration]] <- stage3$safe_auroc(
      heldout_truth[index],
      heldout_probability[index]
    )

    gap[[iteration]] <-
      train_auc[[iteration]] -
      heldout_auc[[iteration]]
  }

  summary <- data.table(
    Development_fold_AUROC_count =
      length(development_fold_aurocs),
    Bootstrap_iterations = iterations,
    Mean_development_AUROC = mean(
      development_fold_aurocs
    ),
    Heldout_AUROC = stage3$safe_auroc(
      heldout_truth,
      heldout_probability
    ),
    Observed_gap =
      mean(development_fold_aurocs) -
      stage3$safe_auroc(
        heldout_truth,
        heldout_probability
      ),
    Gap_lower95 = as.numeric(
      stats::quantile(
        gap,
        0.025,
        names = FALSE
      )
    ),
    Gap_upper95 = as.numeric(
      stats::quantile(
        gap,
        0.975,
        names = FALSE
      )
    ),
    Support_gap_gt_0 = mean(
      gap > 0
    )
  )

  list(
    summary = summary,
    distribution = data.table(
      iteration = seq_len(iterations),
      development_AUROC = train_auc,
      heldout_AUROC = heldout_auc,
      Delta_AUROC = gap
    )
  )
}

run_experiment_1_testing <- function() {
  diseases <- intersect(
    c(
      "TB",
      "RA",
      "UC",
      "CD"
    ),
    REQUESTED_DISEASES
  )

  output_root <- file.path(
    STAGE4_ROOT,
    "01_experiment_1_sample_support"
  )

  make_dir(
    output_root
  )

  all_exact <- list()
  all_delong <- list()
  all_bootstrap <- list()
  counter_exact <- 0L
  counter_delong <- 0L
  counter_bootstrap <- 0L

  for (disease_id in diseases) {
    message_section(
      paste0(
        "Stage 04 | Experiment 1 | ",
        disease_id
      )
    )

    split <- load_equal_size_split(
      disease_id
    )

    disease_dir <- file.path(
      output_root,
      disease_id
    )

    make_dir(
      disease_dir
    )

    safe_write_csv(
      split$audit,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_frozen_split_audit.csv"
        )
      )
    )

    representation_ids <- as.character(
      get_config_value(
        config,
        "sample_support_experiment.representations",
        c(
          "VST_fSVA",
          "VST_ComBatRef",
          "ComBatSeq_log1p"
        )
      )
    )

    prediction_rows <- list()
    metric_rows <- list()

    for (representation_index in seq_along(
      representation_ids
    )) {
      representation_id <- representation_ids[[representation_index]]

      seed <- make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment1",
        disease_id,
        representation_id
      )

      representation <- stage1$fit_apply_representation(
        representation_id = representation_id,
        train_counts = split$development$counts,
        new_counts = split$held_out$counts,
        train_metadata = split$development$metadata,
        new_metadata = split$held_out$metadata,
        seed = seed
      )

      ranking <- fit_rf_importance_ranking(
        train_matrix = representation$train,
        train_metadata = split$development$metadata,
        trees = 300L,
        top_n = 50L,
        seed = seed + 1L
      )

      selected_genes <- ranking$genes

      fit <- fit_fixed_mlr3_model(
        learner_id = "LightGBM",
        train_matrix = representation$train[
          selected_genes,
          ,
          drop = FALSE
        ],
        train_metadata = split$development$metadata,
        seed = seed + 2L
      )

      probability <- predict_fixed_mlr3_model(
        fit,
        representation$new[
          selected_genes,
          ,
          drop = FALSE
        ]
      )

      current_predictions <- data.table(
        disease = disease_id,
        representation = representation_id,
        sample_id = split$held_out$metadata$sample_id,
        truth = as.character(
          split$held_out$metadata$outcome
        ),
        probability = probability
      )

      prediction_rows[[representation_index]] <- current_predictions

      bundle <- save_performance_bundle(
        output_dir = file.path(
          disease_dir,
          representation_id
        ),
        prefix = paste0(
          disease_id,
          "_",
          representation_id,
          "_heldout"
        ),
        truth = split$held_out$metadata$outcome,
        probability = probability,
        threshold = 0.50,
        calibration_bins = 6L,
        bootstrap_iterations = 0L,
        seed = seed,
        save_bootstrap = FALSE
      )

      current_metrics <- data.table::copy(
        bundle$metrics
      )

      current_metrics[
        ,
        `:=`(
          disease = disease_id,
          representation = representation_id
        )
      ]

      metric_rows[[representation_index]] <- current_metrics

      exact <- exact_balanced_permutation_test(
        probability = probability,
        observed_truth = split$held_out$metadata$outcome
      )

      counter_exact <- counter_exact + 1L

      exact$result[
        ,
        `:=`(
          disease = disease_id,
          representation = representation_id
        )
      ]

      all_exact[[counter_exact]] <- exact$result

      safe_write_csv(
        exact$null,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_exact_null_distribution.csv"
          )
        )
      )

      stage3_fold_path <- file.path(
        STAGE3_ROOT,
        "01_experiment_1_sample_support",
        disease_id,
        paste0(
          disease_id,
          "_development_fold_metrics.csv"
        )
      )

      check_file_exists(
        stage3_fold_path,
        paste0(
          disease_id,
          " Stage-03 Experiment-1 fold metrics"
        )
      )

      stage3_folds <- as.data.table(
        safe_read_csv(
          stage3_fold_path
        )
      )

      reference_folds <- stage3_folds[
        representation == representation_id &
          selector == "importance" &
          top_n == 50L &
          learner == "LightGBM"
      ]

      development_fold_aurocs <- reference_folds$AUROC

      fold_audit <- data.table(
        disease = disease_id,
        representation = representation_id,
        methods_text_expected_development_folds = 5L,
        available_stage03_development_folds =
          length(development_fold_aurocs),
        stage03_design =
          "10-fold x 2-repeat CV; all saved matching fold AUROCs retained"
      )

      safe_write_csv(
        fold_audit,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_bootstrap_fold_count_audit.csv"
          )
        )
      )

      if (length(development_fold_aurocs) != 5L) {
        warning(
          disease_id,
          " / ",
          representation_id,
          ": dissertation prose references five development-fold AUROCs, while Stage 03 contains ",
          length(development_fold_aurocs),
          " matching fold results. Stage 04 uses the actual saved Stage-03 values."
        )
      }

      gap_bootstrap <- bootstrap_development_to_heldout_gap(
        development_fold_aurocs = development_fold_aurocs,
        heldout_truth = split$held_out$metadata$outcome,
        heldout_probability = probability,
        iterations = 5000L,
        seed = seed + 3L
      )

      counter_bootstrap <- counter_bootstrap + 1L

      gap_bootstrap$summary[
        ,
        `:=`(
          disease = disease_id,
          representation = representation_id
        )
      ]

      all_bootstrap[[counter_bootstrap]] <- gap_bootstrap$summary

      safe_write_csv(
        gap_bootstrap$summary,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_development_to_heldout_bootstrap_summary.csv"
          )
        )
      )

      safe_write_csv(
        gap_bootstrap$distribution,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_development_to_heldout_bootstrap_distribution.csv"
          )
        )
      )

      safe_save_rds(
        selected_genes,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_final_top50_genes.rds"
          )
        )
      )

      safe_save_rds(
        fit,
        file.path(
          disease_dir,
          representation_id,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_LightGBM_model.rds"
          )
        )
      )
    }

    disease_predictions <- data.table::rbindlist(
      prediction_rows,
      use.names = TRUE,
      fill = TRUE
    )

    disease_metrics <- data.table::rbindlist(
      metric_rows,
      use.names = TRUE,
      fill = TRUE
    )

    safe_write_csv(
      disease_predictions,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_all_representation_heldout_predictions.csv"
        )
      )
    )

    safe_write_csv(
      disease_metrics,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_all_representation_heldout_metrics.csv"
        )
      )
    )

    pairs <- utils::combn(
      representation_ids,
      2L,
      simplify = FALSE
    )

    disease_delong <- vector(
      "list",
      length(pairs)
    )

    for (pair_index in seq_along(pairs)) {
      pair <- pairs[[pair_index]]

      pred_a <- disease_predictions[
        representation == pair[[1]]
      ][
        match(
          split$held_out$metadata$sample_id,
          sample_id
        ),
        probability
      ]

      pred_b <- disease_predictions[
        representation == pair[[2]]
      ][
        match(
          split$held_out$metadata$sample_id,
          sample_id
        ),
        probability
      ]

      current <- paired_delong_test(
        truth = split$held_out$metadata$outcome,
        probability_a = pred_a,
        probability_b = pred_b,
        label_a = pair[[1]],
        label_b = pair[[2]]
      )

      current[
        ,
        disease := disease_id
      ]

      disease_delong[[pair_index]] <- current
    }

    disease_delong <- data.table::rbindlist(
      disease_delong,
      use.names = TRUE,
      fill = TRUE
    )

    disease_delong[
      ,
      p_BH_within_disease := stats::p.adjust(
        p_value,
        method = "BH"
      )
    ]

    safe_write_csv(
      disease_delong,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_paired_DeLong_representation_comparisons.csv"
        )
      )
    )

    for (row_index in seq_len(nrow(disease_delong))) {
      counter_delong <- counter_delong + 1L
      all_delong[[counter_delong]] <- disease_delong[
        row_index
      ]
    }
  }

  exact_dt <- data.table::rbindlist(
    all_exact,
    use.names = TRUE,
    fill = TRUE
  )

  exact_dt[
    ,
    p_BH_global := stats::p.adjust(
      exact_one_sided_p,
      method = "BH"
    )
  ]

  delong_dt <- data.table::rbindlist(
    all_delong,
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(delong_dt) > 0L) {
    delong_dt[
      ,
      p_BH_global := stats::p.adjust(
        p_value,
        method = "BH"
      )
    ]
  }

  bootstrap_dt <- data.table::rbindlist(
    all_bootstrap,
    use.names = TRUE,
    fill = TRUE
  )

  safe_write_csv(
    exact_dt,
    file.path(
      output_root,
      "Experiment1_exact_permutation_tests_BH.csv"
    )
  )

  safe_write_csv(
    delong_dt,
    file.path(
      output_root,
      "Experiment1_paired_DeLong_tests_BH.csv"
    )
  )

  safe_write_csv(
    bootstrap_dt,
    file.path(
      output_root,
      "Experiment1_development_to_heldout_bootstrap_summary.csv"
    )
  )

  invisible(
    list(
      exact = exact_dt,
      delong = delong_dt,
      bootstrap = bootstrap_dt
    )
  )
}



# ==============================================================================
# 8. EXPERIMENT 2 - CROSS-DISEASE TRANSFER AND FROZEN TB TEST
# ==============================================================================

read_stage3_cross_disease_lock <- function(model_name) {
  path <- file.path(
    STAGE3_ROOT,
    "02_experiment_2_disease_composition",
    model_name,
    paste0(
      model_name,
      "_LOCKED_PIPELINE.csv"
    )
  )

  check_file_exists(
    path,
    paste0(
      model_name,
      " Stage-03 cross-disease lock"
    )
  )

  lock <- as.data.table(
    safe_read_csv(
      path
    )
  )

  if (nrow(lock) != 1L) {
    stop(
      model_name,
      " cross-disease lock must contain exactly one row.",
      call. = FALSE
    )
  }

  expected <- c(
    representation = "VST_fSVA",
    selector = "mrmr",
    class_handling = "None",
    learner = "LightGBM"
  )

  audit <- data.table(
    field = names(expected),
    expected = unname(expected),
    observed = c(
      as.character(lock$representation[[1]]),
      as.character(lock$selector[[1]]),
      as.character(lock$class_handling[[1]]),
      as.character(lock$learner[[1]])
    )
  )

  audit[
    ,
    match := expected == observed
  ]

  safe_write_csv(
    audit,
    file.path(
      AUDIT_ROOT,
      paste0(
        "Experiment2_",
        model_name,
        "_lock_audit.csv"
      )
    )
  )

  if (!all(audit$match)) {
    stop(
      model_name,
      " Stage-03 lock differs from the prespecified VST+fSVA | mRMR | None | LightGBM pipeline.",
      call. = FALSE
    )
  }

  lock
}

cross_disease_lock_parameters <- function(lock) {
  required <- c(
    "learning_rate",
    "num_leaves",
    "feature_fraction",
    "bagging_fraction",
    "min_data_in_leaf",
    "num_iterations"
  )

  missing <- setdiff(
    required,
    names(lock)
  )

  if (length(missing) > 0L) {
    stop(
      "Cross-disease LightGBM lock is missing parameter(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  as.list(
    lock[
      1L,
      ..required
    ]
  )
}

fit_tunable_mlr3_model <- function(learner_id,
                                   parameters,
                                   train_matrix,
                                   train_metadata,
                                   seed,
                                   class_weights = NULL) {
  task_info <- stage3$make_mlr3_task(
    expression_matrix = train_matrix,
    metadata = train_metadata,
    task_id = paste0(
      "stage04_tuned_",
      learner_id
    )
  )

  learner <- stage3$make_tunable_learner(
    learner_id = learner_id,
    parameters = parameters,
    seed = seed,
    class_weights = class_weights
  )

  set.seed(
    seed
  )

  learner$train(
    task_info$task
  )

  list(
    learner = learner,
    gene_map = task_info$map
  )
}

fit_cross_disease_locked_system <- function(development_counts,
                                            development_metadata,
                                            evaluation_counts,
                                            evaluation_metadata,
                                            parameters,
                                            seed) {
  representation <- stage1$fit_apply_vst_fsva(
    train_counts = development_counts,
    new_counts = evaluation_counts,
    train_metadata = development_metadata,
    seed = make_trace_subseed(
      seed,
      "VST_fSVA"
    ),
    include_cohort_in_vst_design = FALSE
  )

  selected <- stage2$select_ranked_features(
    train_matrix = representation$train,
    train_metadata = development_metadata,
    filter_id = "mrmr",
    top_n = 50L,
    seed = make_trace_subseed(
      seed,
      "mrmr"
    ),
    rf_trees = 300L
  )

  genes <- selected$genes

  model <- fit_tunable_mlr3_model(
    learner_id = "LightGBM",
    parameters = parameters,
    train_matrix = representation$train[
      genes,
      ,
      drop = FALSE
    ],
    train_metadata = development_metadata,
    seed = make_trace_subseed(
      seed,
      "LightGBM"
    )
  )

  probability <- predict_fixed_mlr3_model(
    model,
    representation$new[
      genes,
      ,
      drop = FALSE
    ]
  )

  list(
    representation = representation,
    selected = selected,
    genes = genes,
    model = model,
    probability = probability
  )
}

run_experiment_2_lodo <- function(common_genes,
                                  parameters,
                                  output_root) {
  disease_ids <- c(
    "TB",
    "RA",
    "UC",
    "CD",
    "T1D"
  )

  prediction_rows <- list()
  metric_rows <- list()

  for (target_index in seq_along(disease_ids)) {
    target_disease <- disease_ids[[target_index]]

    source_diseases <- setdiff(
      disease_ids,
      target_disease
    )

    source <- stage3$prepare_cross_disease_development(
      source_diseases,
      common_genes
    )

    target_canonical <- stage3$load_canonical_disease(
      target_disease
    )

    target_counts <- target_canonical$development_counts[
      common_genes,
      ,
      drop = FALSE
    ]

    target_metadata <- normalise_binary_metadata(
      target_canonical$development_metadata
    )

    seed <- make_trace_subseed(
      GLOBAL_SEED,
      "stage04",
      "experiment2",
      "LODO",
      target_disease
    )

    fitted <- fit_cross_disease_locked_system(
      development_counts = source$counts,
      development_metadata = source$metadata,
      evaluation_counts = target_counts,
      evaluation_metadata = target_metadata,
      parameters = parameters,
      seed = seed
    )

    current_predictions <- data.table(
      target_disease = target_disease,
      source_diseases = paste(
        source_diseases,
        collapse = "+"
      ),
      sample_id = target_metadata$sample_id,
      truth = as.character(
        target_metadata$outcome
      ),
      probability = fitted$probability
    )

    current_metrics <- performance_metrics(
      target_metadata$outcome,
      fitted$probability,
      threshold = 0.50
    )

    current_metrics[
      ,
      `:=`(
        target_disease = target_disease,
        source_diseases = paste(
          source_diseases,
          collapse = "+"
        ),
        target_n = nrow(target_metadata)
      )
    ]

    prediction_rows[[target_index]] <- current_predictions
    metric_rows[[target_index]] <- current_metrics

    target_dir <- file.path(
      output_root,
      "LODO",
      target_disease
    )

    make_dir(
      target_dir
    )

    safe_write_csv(
      current_predictions,
      file.path(
        target_dir,
        paste0(
          "LODO_",
          target_disease,
          "_predictions.csv"
        )
      )
    )

    safe_write_csv(
      current_metrics,
      file.path(
        target_dir,
        paste0(
          "LODO_",
          target_disease,
          "_metrics.csv"
        )
      )
    )

    safe_save_rds(
      fitted$genes,
      file.path(
        target_dir,
        paste0(
          "LODO_",
          target_disease,
          "_mRMR_top50_genes.rds"
        )
      )
    )
  }

  predictions <- data.table::rbindlist(
    prediction_rows,
    use.names = TRUE,
    fill = TRUE
  )

  metrics <- data.table::rbindlist(
    metric_rows,
    use.names = TRUE,
    fill = TRUE
  )

  safe_write_csv(
    predictions,
    file.path(
      output_root,
      "LODO_all_target_predictions.csv"
    )
  )

  safe_write_csv(
    metrics,
    file.path(
      output_root,
      "LODO_all_target_metrics.csv"
    )
  )

  list(
    predictions = predictions,
    metrics = metrics
  )
}

run_experiment_2_testing <- function() {
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

  message_section(
    "Stage 04 | Experiment 2 | cross-disease transfer and frozen TB testing"
  )

  output_root <- file.path(
    STAGE4_ROOT,
    "02_experiment_2_disease_composition"
  )

  make_dir(
    output_root
  )

  common_genes <- stage3$load_cross_disease_common_genes()

  if (length(common_genes) != 4450L) {
    warning(
      "Observed cross-disease common gene count is ",
      length(common_genes),
      "; dissertation reference = 4,450."
    )

    if (STRICT_AUDITS) {
      stop(
        "Strict audit requested and common cross-disease gene count differs from 4,450.",
        call. = FALSE
      )
    }
  }

  model_a_lock <- read_stage3_cross_disease_lock(
    "Model_A"
  )

  model_b_lock <- read_stage3_cross_disease_lock(
    "Model_B"
  )

  model_a_parameters <- cross_disease_lock_parameters(
    model_a_lock
  )

  model_b_parameters <- cross_disease_lock_parameters(
    model_b_lock
  )

  # The locked LODO analytical architecture is fixed. The classifier
  # configuration below is sourced from the all-five development lock because
  # Stage 03 stores final model-specific hyperparameters after development-only
  # tuning. This parameter source is explicitly audited.
  lodo_parameter_audit <- data.table(
    analytical_pipeline =
      "VST_fSVA | mRMR | top50 | None | LightGBM",
    LightGBM_parameter_source =
      "Stage03 Model_B all-five development-only lock",
    final_target_labels_used = FALSE
  )

  safe_write_csv(
    lodo_parameter_audit,
    file.path(
      output_root,
      "LODO_parameter_source_audit.csv"
    )
  )

  lodo <- run_experiment_2_lodo(
    common_genes = common_genes,
    parameters = model_b_parameters,
    output_root = output_root
  )

  model_a <- stage3$prepare_cross_disease_development(
    c(
      "RA",
      "UC",
      "CD",
      "T1D"
    ),
    common_genes
  )

  model_b <- stage3$prepare_cross_disease_development(
    c(
      "TB",
      "RA",
      "UC",
      "CD",
      "T1D"
    ),
    common_genes
  )

  tb_split <- load_general_split(
    "TB"
  )

  frozen_tb_counts <- tb_split$held_out$counts[
    common_genes,
    ,
    drop = FALSE
  ]

  frozen_tb_metadata <- tb_split$held_out$metadata

  audit_expected_cross_disease_counts(
    model_a_n = nrow(model_a$metadata),
    model_b_n = nrow(model_b$metadata),
    tb_test_n = nrow(frozen_tb_metadata)
  )

  final_results <- list()
  final_predictions <- list()
  final_importance <- list()

  model_specifications <- list(
    Model_A = list(
      data = model_a,
      lock = model_a_lock,
      parameters = model_a_parameters,
      label = "Model A: TB-naive"
    ),
    Model_B = list(
      data = model_b,
      lock = model_b_lock,
      parameters = model_b_parameters,
      label = "Model B: ALL5"
    )
  )

  for (model_index in seq_along(
    model_specifications
  )) {
    model_name <- names(
      model_specifications
    )[[model_index]]

    specification <- model_specifications[[model_name]]

    model_dir <- file.path(
      output_root,
      "frozen_TB",
      model_name
    )

    make_dir(
      model_dir
    )

    seed <- make_trace_subseed(
      GLOBAL_SEED,
      "stage04",
      "experiment2",
      "frozen_TB",
      model_name
    )

    fitted <- fit_cross_disease_locked_system(
      development_counts = specification$data$counts,
      development_metadata = specification$data$metadata,
      evaluation_counts = frozen_tb_counts,
      evaluation_metadata = frozen_tb_metadata,
      parameters = specification$parameters,
      seed = seed
    )

    threshold <- as.numeric(
      specification$lock$development_threshold[[1L]]
    )

    bundle <- save_performance_bundle(
      output_dir = model_dir,
      prefix = paste0(
        model_name,
        "_frozen_TB"
      ),
      truth = frozen_tb_metadata$outcome,
      probability = fitted$probability,
      threshold = threshold,
      calibration_bins = 6L,
      bootstrap_iterations = 2000L,
      seed = seed,
      save_bootstrap = TRUE
    )

    predictions <- data.table(
      Model = model_name,
      Model_label = specification$label,
      sample_id = frozen_tb_metadata$sample_id,
      cohort_group = as.character(
        frozen_tb_metadata$cohort_group
      ),
      truth = as.character(
        frozen_tb_metadata$outcome
      ),
      probability = fitted$probability,
      threshold = threshold
    )

    safe_write_csv(
      predictions,
      file.path(
        model_dir,
        paste0(
          model_name,
          "_frozen_TB_predictions.csv"
        )
      )
    )

    safe_save_rds(
      fitted$model,
      file.path(
        model_dir,
        paste0(
          model_name,
          "_final_LightGBM_model.rds"
        )
      )
    )

    safe_save_rds(
      fitted$genes,
      file.path(
        model_dir,
        paste0(
          model_name,
          "_final_mRMR_top50_genes.rds"
        )
      )
    )

    safe_save_rds(
      list(
        representation_state = fitted$representation$state,
        genes = fitted$genes,
        threshold = threshold,
        LightGBM_parameters = specification$parameters,
        development_sample_ids =
          specification$data$metadata$sample_id,
        frozen_TB_sample_ids =
          frozen_tb_metadata$sample_id
      ),
      file.path(
        model_dir,
        paste0(
          model_name,
          "_final_state.rds"
        )
      )
    )

    x_heldout <- as.data.frame(
      t(
        fitted$representation$new[
          fitted$genes,
          ,
          drop = FALSE
        ]
      ),
      check.names = FALSE
    )

    predict_probability <- function(newdata) {
      matrix_new <- t(
        as.matrix(
          newdata
        )
      )

      rownames(matrix_new) <- colnames(
        newdata
      )

      colnames(matrix_new) <- rownames(
        newdata
      )

      predict_fixed_mlr3_model(
        fitted$model,
        matrix_new
      )
    }

    importance <- run_delta_auroc_permutation(
      model_label = model_name,
      x = x_heldout,
      truth = frozen_tb_metadata$outcome,
      predict_probability = predict_probability,
      permutations = 50L,
      seed = seed + 10L,
      output_dir = file.path(
        model_dir,
        "permutation_importance"
      ),
      resume = TRUE
    )

    reactome <- run_reactome_enrichment(
      importance_summary = importance$summary,
      background_genes = common_genes,
      output_dir = file.path(
        model_dir,
        "reactome"
      ),
      prefix = paste0(
        model_name,
        "_frozen_TB"
      ),
      top_positive_genes = 25L
    )

    final_predictions[[model_name]] <- predictions

    current_metrics <- data.table::copy(
      bundle$metrics
    )

    current_metrics[
      ,
      `:=`(
        Model = model_name,
        Model_label = specification$label
      )
    ]

    final_results[[model_name]] <- current_metrics
    final_importance[[model_name]] <- importance$summary
  }

  model_a_prediction <- final_predictions$Model_A[
    match(
      frozen_tb_metadata$sample_id,
      sample_id
    ),
    probability
  ]

  model_b_prediction <- final_predictions$Model_B[
    match(
      frozen_tb_metadata$sample_id,
      sample_id
    ),
    probability
  ]

  comparison <- paired_delong_test(
    truth = frozen_tb_metadata$outcome,
    probability_a = model_a_prediction,
    probability_b = model_b_prediction,
    label_a = "Model A: TB-naive",
    label_b = "Model B: ALL5"
  )

  comparison[
    ,
    Delta_ModelB_minus_ModelA :=
      AUROC_B -
        AUROC_A
  ]

  safe_write_csv(
    comparison,
    file.path(
      output_root,
      "Frozen_TB_ModelA_vs_ModelB_paired_DeLong.csv"
    )
  )

  metrics <- data.table::rbindlist(
    final_results,
    use.names = TRUE,
    fill = TRUE
  )

  predictions <- data.table::rbindlist(
    final_predictions,
    use.names = TRUE,
    fill = TRUE
  )

  safe_write_csv(
    metrics,
    file.path(
      output_root,
      "Frozen_TB_ModelA_ModelB_metrics.csv"
    )
  )

  safe_write_csv(
    predictions,
    file.path(
      output_root,
      "Frozen_TB_ModelA_ModelB_predictions.csv"
    )
  )

  importance_overlap <- merge(
    final_importance$Model_A[
      ,
      .(
        Gene,
        Model_A_Delta_AUROC = Mean_Delta_AUROC,
        Model_A_Rank = Rank
      )
    ],
    final_importance$Model_B[
      ,
      .(
        Gene,
        Model_B_Delta_AUROC = Mean_Delta_AUROC,
        Model_B_Rank = Rank
      )
    ],
    by = "Gene",
    all = TRUE
  )

  safe_write_csv(
    importance_overlap,
    file.path(
      output_root,
      "Frozen_TB_ModelA_ModelB_gene_importance_alignment.csv"
    )
  )

  invisible(
    list(
      lodo = lodo,
      final_metrics = metrics,
      final_predictions = predictions,
      paired_delong = comparison,
      importance = final_importance
    )
  )
}



# ==============================================================================
# 9. EXPERIMENT 3 - COHORT TRANSPORTABILITY, FINAL TEST AND RELIABILITY
# ==============================================================================

read_stage3_disease_lock <- function(disease_id) {
  path <- file.path(
    STAGE3_ROOT,
    "03_experiment_3_cohort_composition",
    disease_id,
    paste0(
      disease_id,
      "_LOCKED_DISEASE_SPECIFIC_PIPELINE.csv"
    )
  )

  check_file_exists(
    path,
    paste0(
      disease_id,
      " Stage-03 disease-specific lock"
    )
  )

  lock <- as.data.table(
    safe_read_csv(
      path
    )
  )

  if (nrow(lock) != 1L) {
    stop(
      disease_id,
      " disease-specific Stage-03 lock must contain exactly one row.",
      call. = FALSE
    )
  }

  required <- c(
    "representation",
    "weighting",
    "top_n",
    "model",
    "development_threshold"
  )

  missing <- setdiff(
    required,
    names(lock)
  )

  if (length(missing) > 0L) {
    stop(
      disease_id,
      " lock is missing column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (lock$model[[1L]] != "Random_Forest") {
    stop(
      disease_id,
      " final Stage-04 disease-specific implementation expects the Stage-03 Random Forest lock.",
      call. = FALSE
    )
  }

  lock
}

load_stage3_selected_loco_predictions <- function(disease_id,
                                                  lock) {
  path <- file.path(
    STAGE3_ROOT,
    "03_experiment_3_cohort_composition",
    disease_id,
    paste0(
      disease_id,
      "_development_LOCO_predictions_all_candidates.csv"
    )
  )

  check_file_exists(
    path,
    paste0(
      disease_id,
      " Stage-03 LOCO predictions"
    )
  )

  predictions <- as.data.table(
    safe_read_csv(
      path
    )
  )

  selected <- predictions[
    representation ==
      as.character(lock$representation[[1L]]) &
      weighting ==
        as.character(lock$weighting[[1L]]) &
      as.integer(top_n) ==
        as.integer(lock$top_n[[1L]]) &
      model ==
        "Random_Forest"
  ]

  if (nrow(selected) == 0L) {
    stop(
      disease_id,
      " selected Stage-03 LOCO configuration has no prediction rows.",
      call. = FALSE
    )
  }

  selected
}

summarise_selected_loco <- function(disease_id,
                                    selected_predictions,
                                    threshold,
                                    output_dir) {
  target_cohorts <- unique(
    as.character(
      selected_predictions$target_cohort
    )
  )

  metric_rows <- list()

  for (cohort_index in seq_along(target_cohorts)) {
    target_cohort <- target_cohorts[[cohort_index]]

    current <- selected_predictions[
      get("target_cohort") ==
        target_cohort
    ]

    current_truth <- factor(
      current$truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    )

    cohort_dir <- file.path(
      output_dir,
      "LOCO",
      safe_file_name(
        target_cohort
      )
    )

    bundle <- save_performance_bundle(
      output_dir = cohort_dir,
      prefix = paste0(
        disease_id,
        "_LOCO_",
        safe_file_name(target_cohort)
      ),
      truth = current_truth,
      probability = current$probability,
      threshold = threshold,
      calibration_bins = 6L,
      bootstrap_iterations = 2000L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment3",
        disease_id,
        "LOCO",
        target_cohort
      ),
      save_bootstrap = TRUE
    )

    current_metrics <- data.table::copy(
      bundle$metrics
    )

    current_metrics[
      ,
      `:=`(
        disease = disease_id,
        target_cohort = target_cohort,
        target_n = nrow(current)
      )
    ]

    metric_rows[[cohort_index]] <- current_metrics

    safe_write_csv(
      current,
      file.path(
        cohort_dir,
        paste0(
          disease_id,
          "_",
          safe_file_name(target_cohort),
          "_selected_LOCO_predictions.csv"
        )
      )
    )
  }

  metrics <- data.table::rbindlist(
    metric_rows,
    use.names = TRUE,
    fill = TRUE
  )

  summary <- metrics[
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
      )
    )
  ]

  safe_write_csv(
    metrics,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_selected_pipeline_LOCO_metrics.csv"
      )
    )
  )

  safe_write_csv(
    summary,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_selected_pipeline_LOCO_summary.csv"
      )
    )
  )

  list(
    metrics = metrics,
    summary = summary
  )
}

fit_final_disease_specific_system <- function(disease_id,
                                              split,
                                              lock,
                                              output_dir) {
  representation_id <- as.character(
    lock$representation[[1L]]
  )

  top_n <- as.integer(
    lock$top_n[[1L]]
  )

  weighting <- as.character(
    lock$weighting[[1L]]
  )

  threshold <- as.numeric(
    lock$development_threshold[[1L]]
  )

  seed <- make_trace_subseed(
    GLOBAL_SEED,
    "stage04",
    "experiment3",
    disease_id,
    "final_heldout"
  )

  representation <- stage1$fit_apply_representation(
    representation_id = representation_id,
    train_counts = split$development$counts,
    new_counts = split$held_out$counts,
    train_metadata = split$development$metadata,
    new_metadata = split$held_out$metadata,
    seed = seed
  )

  ranking <- fit_rf_importance_ranking(
    train_matrix = representation$train,
    train_metadata = split$development$metadata,
    trees = 1000L,
    top_n = top_n,
    seed = seed + 1L
  )

  genes <- ranking$genes

  train_matrix <- representation$train[
    genes,
    ,
    drop = FALSE
  ]

  heldout_matrix <- representation$new[
    genes,
    ,
    drop = FALSE
  ]

  model <- fit_final_ranger_model(
    train_matrix = train_matrix,
    train_metadata = split$development$metadata,
    weighting = weighting,
    seed = seed + 2L,
    trees = 250L,
    min_node_size = 5L
  )

  probability <- predict_final_ranger_model(
    model,
    heldout_matrix
  )

  bundle <- save_performance_bundle(
    output_dir = file.path(
      output_dir,
      "frozen_heldout"
    ),
    prefix = paste0(
      disease_id,
      "_locked_RF_heldout"
    ),
    truth = split$held_out$metadata$outcome,
    probability = probability,
    threshold = threshold,
    calibration_bins = 6L,
    bootstrap_iterations = 2000L,
    seed = seed + 3L,
    save_bootstrap = TRUE
  )

  predictions <- data.table(
    disease = disease_id,
    sample_id = split$held_out$metadata$sample_id,
    cohort_group = as.character(
      split$held_out$metadata$cohort_group
    ),
    truth = as.character(
      split$held_out$metadata$outcome
    ),
    probability = probability,
    threshold = threshold,
    representation = representation_id,
    top_n = top_n,
    weighting = weighting,
    model = "Random_Forest"
  )

  safe_write_csv(
    predictions,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_locked_RF_heldout_predictions.csv"
      )
    )
  )

  safe_save_rds(
    ranking$scores,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_final_RF_importance_ranking.rds"
      )
    )
  )

  safe_save_rds(
    genes,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_final_selected_genes.rds"
      )
    )
  )

  model_path <- file.path(
    output_dir,
    "frozen_heldout",
    paste0(
      disease_id,
      "_final_RF_model.rds"
    )
  )

  safe_save_rds(
    model,
    model_path
  )

  safe_save_rds(
    train_matrix,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_final_development_matrix.rds"
      )
    )
  )

  safe_save_rds(
    heldout_matrix,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_final_heldout_matrix.rds"
      )
    )
  )

  restored_model <- safe_read_rds(
    model_path
  )

  restored_probability <- predict_final_ranger_model(
    restored_model,
    heldout_matrix
  )

  max_difference <- max(
    abs(
      restored_probability -
        probability
    ),
    na.rm = TRUE
  )

  restoration_audit <- data.table(
    disease = disease_id,
    model = "Random_Forest",
    heldout_n = length(probability),
    max_absolute_probability_difference =
      max_difference,
    tolerance = 1e-12,
    passed = max_difference <= 1e-12
  )

  safe_write_csv(
    restoration_audit,
    file.path(
      output_dir,
      "frozen_heldout",
      paste0(
        disease_id,
        "_RF_restoration_audit.csv"
      )
    )
  )

  if (!restoration_audit$passed[[1]]) {
    stop(
      disease_id,
      " restored Random Forest predictions differ from the original probabilities.",
      call. = FALSE
    )
  }

  list(
    representation = representation,
    ranking = ranking,
    genes = genes,
    model = model,
    train_matrix = train_matrix,
    heldout_matrix = heldout_matrix,
    probability = probability,
    predictions = predictions,
    bundle = bundle,
    restoration_audit = restoration_audit
  )
}

run_experiment_3_testing <- function() {
  output_root <- file.path(
    STAGE4_ROOT,
    "03_experiment_3_cohort_composition"
  )

  make_dir(
    output_root
  )

  overall_rows <- list()
  loco_rows <- list()
  output_counter <- 0L

  for (disease_id in REQUESTED_DISEASES) {
    message_section(
      paste0(
        "Stage 04 | Experiment 3 | ",
        disease_id
      )
    )

    disease_dir <- file.path(
      output_root,
      disease_id
    )

    make_dir(
      disease_dir
    )

    lock <- read_stage3_disease_lock(
      disease_id
    )

    safe_write_csv(
      lock,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_Stage03_lock_used_by_Stage04.csv"
        )
      )
    )

    selected_loco <- load_stage3_selected_loco_predictions(
      disease_id,
      lock
    )

    threshold <- as.numeric(
      lock$development_threshold[[1L]]
    )

    loco <- summarise_selected_loco(
      disease_id = disease_id,
      selected_predictions = selected_loco,
      threshold = threshold,
      output_dir = disease_dir
    )

    split <- load_general_split(
      disease_id
    )

    safe_write_csv(
      split$audit,
      file.path(
        disease_dir,
        paste0(
          disease_id,
          "_general_partition_audit.csv"
        )
      )
    )

    final <- fit_final_disease_specific_system(
      disease_id = disease_id,
      split = split,
      lock = lock,
      output_dir = disease_dir
    )

    output_counter <- output_counter + 1L

    current_final <- data.table::copy(
      final$bundle$metrics
    )

    current_final[
      ,
      `:=`(
        disease = disease_id,
        representation =
          as.character(lock$representation[[1L]]),
        top_n =
          as.integer(lock$top_n[[1L]]),
        weighting =
          as.character(lock$weighting[[1L]])
      )
    ]

    overall_rows[[output_counter]] <- current_final

    current_loco <- data.table::copy(
      loco$summary
    )

    current_loco[
      ,
      disease := disease_id
    ]

    loco_rows[[output_counter]] <- current_loco

    x_heldout <- as.data.frame(
      t(
        final$heldout_matrix
      ),
      check.names = FALSE
    )

    ranger_predict_probability <- function(newdata) {
      x <- as.data.frame(
        newdata,
        check.names = FALSE
      )

      colnames(x) <- final$model$safe_genes

      prediction <- predict(
        final$model$model,
        data = x
      )$predictions

      as.numeric(
        prediction[
          ,
          POSITIVE_CLASS
        ]
      )
    }

    importance <- run_delta_auroc_permutation(
      model_label = paste0(
        disease_id,
        "_Random_Forest"
      ),
      x = x_heldout,
      truth = split$held_out$metadata$outcome,
      predict_probability =
        ranger_predict_probability,
      permutations = 100L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment3",
        disease_id,
        "importance"
      ),
      output_dir = file.path(
        disease_dir,
        "molecular_reliability",
        "permutation_importance"
      ),
      resume = TRUE
    )

    top_genes <- importance$summary[
      order(-Mean_Delta_AUROC)
    ][
      seq_len(
        min(
          6L,
          .N
        )
      ),
      Gene
    ]

    safe_write_csv(
      data.table(
        Gene = top_genes,
        Rank = seq_along(
          top_genes
        )
      ),
      file.path(
        disease_dir,
        "molecular_reliability",
        paste0(
          disease_id,
          "_ALE_selected_genes.csv"
        )
      )
    )

    if (
      requireNamespace(
        "DALEX",
        quietly = TRUE
      ) &&
        requireNamespace(
          "ingredients",
          quietly = TRUE
        )
    ) {
      predict_function <- make_ranger_prediction_function(
        final$model
      )

      explainer <- make_dalex_explainer(
        model = final$model$model,
        x = x_heldout,
        truth = split$held_out$metadata$outcome,
        predict_function = predict_function,
        label = paste0(
          disease_id,
          "_Random_Forest"
        )
      )

      run_ale_profiles(
        explainer = explainer,
        genes = top_genes,
        output_dir = file.path(
          disease_dir,
          "molecular_reliability",
          "ALE"
        ),
        prefix = paste0(
          disease_id,
          "_Random_Forest"
        )
      )
    } else {
      warning(
        disease_id,
        ": DALEX/ingredients unavailable; ALE skipped."
      )
    }
  }

  final_metrics <- data.table::rbindlist(
    overall_rows,
    use.names = TRUE,
    fill = TRUE
  )

  loco_summary <- data.table::rbindlist(
    loco_rows,
    use.names = TRUE,
    fill = TRUE
  )

  combined <- merge(
    final_metrics,
    loco_summary,
    by = "disease",
    all = TRUE
  )

  combined[
    ,
    Frozen_minus_mean_LOCO_AUROC :=
      AUROC -
        Mean_LOCO_AUROC
  ]

  combined[
    ,
    Frozen_minus_worst_LOCO_AUROC :=
      AUROC -
        Worst_cohort_AUROC
  ]

  safe_write_csv(
    final_metrics,
    file.path(
      output_root,
      "Experiment3_all_disease_frozen_heldout_metrics.csv"
    )
  )

  safe_write_csv(
    loco_summary,
    file.path(
      output_root,
      "Experiment3_all_disease_LOCO_summary.csv"
    )
  )

  safe_write_csv(
    combined,
    file.path(
      output_root,
      "Experiment3_sample_level_vs_cohort_level_generalisation.csv"
    )
  )

  invisible(
    list(
      frozen = final_metrics,
      loco = loco_summary,
      comparison = combined
    )
  )
}



# ==============================================================================
# 10. EXPERIMENT 4 - ARCHITECTURE TRANSPORTABILITY, FINAL TEST AND RELIABILITY
# ==============================================================================

read_stage3_architecture_lock <- function(disease_id) {
  directory <- file.path(
    STAGE3_ROOT,
    "04_experiment_4_model_architecture",
    disease_id
  )

  lock_path <- file.path(
    directory,
    paste0(
      disease_id,
      "_LOCKED_ARCHITECTURE_DEVELOPMENT_SPECIFICATION.csv"
    )
  )

  genes_path <- file.path(
    directory,
    paste0(
      disease_id,
      "_LOCKED_ARCHITECTURE_TOP50_GENES.rds"
    )
  )

  check_file_exists(
    lock_path,
    paste0(
      disease_id,
      " Stage-03 architecture lock"
    )
  )

  check_file_exists(
    genes_path,
    paste0(
      disease_id,
      " Stage-03 architecture gene lock"
    )
  )

  lock <- as.data.table(
    safe_read_csv(
      lock_path
    )
  )

  genes <- as.character(
    safe_read_rds(
      genes_path
    )
  )

  if (
    nrow(lock) != 1L ||
      length(genes) != 50L
  ) {
    stop(
      disease_id,
      " Stage-03 architecture lock is incomplete.",
      call. = FALSE
    )
  }

  list(
    directory = directory,
    lock = lock,
    genes = genes
  )
}

architecture_lock_parameters <- function(lock) {
  required <- c(
    "n_estimators",
    "learning_rate",
    "num_leaves",
    "max_depth",
    "min_child_samples",
    "subsample",
    "colsample_bytree",
    "reg_alpha",
    "reg_lambda"
  )

  missing <- setdiff(
    required,
    names(lock)
  )

  if (length(missing) > 0L) {
    stop(
      "Architecture lock is missing LightGBM field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  lock[
    ,
    ..required
  ]
}

architecture_source_only_tuning <- function(disease_id,
                                            source_counts,
                                            source_metadata,
                                            target_cohort,
                                            python) {
  inner_folds <- 3L

  fold_id <- stage3$make_joint_stratified_folds(
    metadata = source_metadata,
    number_of_folds = inner_folds,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "stage04",
      "experiment4",
      disease_id,
      target_cohort,
      "LOCO_inner"
    )
  )

  rows <- list()
  counter <- 0L

  for (inner_fold in seq_len(inner_folds)) {
    train_metadata <- source_metadata[
      fold_id != inner_fold,
      ,
      drop = FALSE
    ]

    valid_metadata <- source_metadata[
      fold_id == inner_fold,
      ,
      drop = FALSE
    ]

    train_counts <- source_counts[
      ,
      train_metadata$sample_id,
      drop = FALSE
    ]

    valid_counts <- source_counts[
      ,
      valid_metadata$sample_id,
      drop = FALSE
    ]

    representation <- stage1$fit_apply_log2_cpm_plus_1(
      train_counts = train_counts,
      new_counts = valid_counts
    )

    ranking <- fit_rf_importance_ranking(
      train_matrix = representation$train,
      train_metadata = train_metadata,
      trees = 250L,
      top_n = 50L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment4",
        disease_id,
        target_cohort,
        "LOCO_inner",
        inner_fold,
        "genes"
      )
    )

    genes <- ranking$genes

    x_train <- t(
      representation$train[
        genes,
        ,
        drop = FALSE
      ]
    )

    x_valid <- t(
      representation$new[
        genes,
        ,
        drop = FALSE
      ]
    )

    for (config_index in seq_len(
      nrow(stage3$LIGHTGBM_GRID)
    )) {
      parameters <- stage3$LIGHTGBM_GRID[
        config_index
      ]

      probability <- stage3$fit_predict_python_lightgbm(
        x_train = x_train,
        y_train = train_metadata$outcome,
        x_test = x_valid,
        parameters = parameters,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          "stage04",
          "experiment4",
          disease_id,
          target_cohort,
          "LOCO_inner",
          inner_fold,
          parameters$Config_ID[[1L]]
        ),
        python = python
      )

      counter <- counter + 1L

      rows[[counter]] <- data.table(
        disease = disease_id,
        target_cohort = target_cohort,
        inner_fold = inner_fold,
        Config_ID = parameters$Config_ID[[1L]],
        AUROC = stage3$safe_auroc(
          valid_metadata$outcome,
          probability
        )
      )
    }
  }

  raw <- data.table::rbindlist(
    rows,
    use.names = TRUE,
    fill = TRUE
  )

  summary <- raw[
    ,
    .(
      Inner_folds = .N,
      Mean_inner_AUROC = mean(
        AUROC,
        na.rm = TRUE
      ),
      SD_inner_AUROC = stats::sd(
        AUROC,
        na.rm = TRUE
      ),
      Worst_inner_AUROC = min(
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
    -Worst_inner_AUROC,
    Config_ID
  )

  list(
    best_config_id = summary$Config_ID[[1L]],
    raw = raw,
    summary = summary
  )
}

run_experiment_4_loco <- function(disease_id,
                                  canonical,
                                  python,
                                  output_dir) {
  development_metadata <- normalise_binary_metadata(
    canonical$development_metadata
  )

  target_cohorts <- unique(
    as.character(
      development_metadata$cohort_group
    )
  )

  prediction_rows <- list()
  metric_rows <- list()

  for (target_index in seq_along(target_cohorts)) {
    target_cohort <- target_cohorts[[target_index]]

    source_metadata <- development_metadata[
      as.character(cohort_group) !=
        target_cohort,
      ,
      drop = FALSE
    ]

    target_metadata <- development_metadata[
      as.character(cohort_group) ==
        target_cohort,
      ,
      drop = FALSE
    ]

    source_counts <- canonical$development_counts[
      ,
      source_metadata$sample_id,
      drop = FALSE
    ]

    target_counts <- canonical$development_counts[
      ,
      target_metadata$sample_id,
      drop = FALSE
    ]

    tuning <- architecture_source_only_tuning(
      disease_id = disease_id,
      source_counts = source_counts,
      source_metadata = source_metadata,
      target_cohort = target_cohort,
      python = python
    )

    selected_config <- stage3$LIGHTGBM_GRID[
      Config_ID ==
        tuning$best_config_id
    ]

    representation <- stage1$fit_apply_log2_cpm_plus_1(
      train_counts = source_counts,
      new_counts = target_counts
    )

    ranking <- fit_rf_importance_ranking(
      train_matrix = representation$train,
      train_metadata = source_metadata,
      trees = 250L,
      top_n = 50L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment4",
        disease_id,
        target_cohort,
        "LOCO_full_genes"
      )
    )

    genes <- ranking$genes

    x_train <- t(
      representation$train[
        genes,
        ,
        drop = FALSE
      ]
    )

    x_test <- t(
      representation$new[
        genes,
        ,
        drop = FALSE
      ]
    )

    lgbm_probability <- stage3$fit_predict_python_lightgbm(
      x_train = x_train,
      y_train = source_metadata$outcome,
      x_test = x_test,
      parameters = selected_config,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment4",
        disease_id,
        target_cohort,
        "LOCO_LightGBM"
      ),
      python = python
    )

    tabpfn_probability <- stage3$fit_predict_python_tabpfn(
      x_train = x_train,
      y_train = source_metadata$outcome,
      x_test = x_test,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment4",
        disease_id,
        target_cohort,
        "LOCO_TabPFN"
      ),
      python = python
    )

    target_dir <- file.path(
      output_dir,
      "LOCO",
      safe_file_name(
        target_cohort
      )
    )

    make_dir(
      target_dir
    )

    safe_write_csv(
      tuning$raw,
      file.path(
        target_dir,
        paste0(
          disease_id,
          "_",
          safe_file_name(target_cohort),
          "_LightGBM_source_only_inner_tuning_raw.csv"
        )
      )
    )

    safe_write_csv(
      tuning$summary,
      file.path(
        target_dir,
        paste0(
          disease_id,
          "_",
          safe_file_name(target_cohort),
          "_LightGBM_source_only_inner_tuning_summary.csv"
        )
      )
    )

    safe_save_rds(
      genes,
      file.path(
        target_dir,
        paste0(
          disease_id,
          "_",
          safe_file_name(target_cohort),
          "_source_only_top50_genes.rds"
        )
      )
    )

    for (model_name in c(
      "LightGBM",
      "TabPFN"
    )) {
      probability <- if (
        model_name ==
          "LightGBM"
      ) {
        lgbm_probability
      } else {
        tabpfn_probability
      }

      bundle <- save_performance_bundle(
        output_dir = file.path(
          target_dir,
          model_name
        ),
        prefix = paste0(
          disease_id,
          "_",
          safe_file_name(target_cohort),
          "_",
          model_name,
          "_LOCO"
        ),
        truth = target_metadata$outcome,
        probability = probability,
        threshold = 0.50,
        calibration_bins = 8L,
        bootstrap_iterations = 0L,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          "stage04",
          "experiment4",
          disease_id,
          target_cohort,
          model_name,
          "curves"
        ),
        save_bootstrap = FALSE
      )

      current_predictions <- data.table(
        disease = disease_id,
        target_cohort = target_cohort,
        model = model_name,
        selected_LightGBM_config =
          tuning$best_config_id,
        sample_id = target_metadata$sample_id,
        truth = as.character(
          target_metadata$outcome
        ),
        probability = probability
      )

      current_metrics <- data.table::copy(
        bundle$metrics
      )

      current_metrics[
        ,
        `:=`(
          disease = disease_id,
          target_cohort = target_cohort,
          model = model_name,
          selected_LightGBM_config =
            tuning$best_config_id
        )
      ]

      prediction_rows[[
        length(prediction_rows) + 1L
      ]] <- current_predictions

      metric_rows[[
        length(metric_rows) + 1L
      ]] <- current_metrics
    }
  }

  predictions <- data.table::rbindlist(
    prediction_rows,
    use.names = TRUE,
    fill = TRUE
  )

  metrics <- data.table::rbindlist(
    metric_rows,
    use.names = TRUE,
    fill = TRUE
  )

  summary <- metrics[
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
    predictions,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_architecture_LOCO_predictions.csv"
      )
    )
  )

  safe_write_csv(
    metrics,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_architecture_LOCO_metrics.csv"
      )
    )
  )

  safe_write_csv(
    summary,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_architecture_LOCO_summary.csv"
      )
    )
  )

  list(
    predictions = predictions,
    metrics = metrics,
    summary = summary
  )
}

run_experiment_4_final_disease <- function(disease_id,
                                           python,
                                           output_root) {
  lock_object <- read_stage3_architecture_lock(
    disease_id
  )

  lock <- lock_object$lock
  genes <- lock_object$genes

  split <- load_general_split(
    disease_id
  )

  disease_dir <- file.path(
    output_root,
    disease_id
  )

  make_dir(
    disease_dir
  )

  safe_write_csv(
    split$audit,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_general_partition_audit.csv"
      )
    )
  )

  safe_write_csv(
    lock,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_Stage03_architecture_lock_used_by_Stage04.csv"
      )
    )
  )

  safe_save_rds(
    genes,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_Stage03_locked_top50_genes.rds"
      )
    )
  )

  # Source-only LOCO is development-data transportability evidence.
  loco <- run_experiment_4_loco(
    disease_id = disease_id,
    canonical = split$canonical,
    python = python,
    output_dir = disease_dir
  )

  representation <- stage1$fit_apply_log2_cpm_plus_1(
    train_counts = split$development$counts,
    new_counts = split$held_out$counts
  )

  missing_train <- setdiff(
    genes,
    rownames(representation$train)
  )

  missing_test <- setdiff(
    genes,
    rownames(representation$new)
  )

  if (
    length(missing_train) > 0L ||
      length(missing_test) > 0L
  ) {
    stop(
      disease_id,
      " locked architecture genes are absent from the reconstructed matrix.",
      call. = FALSE
    )
  }

  x_train <- t(
    representation$train[
      genes,
      ,
      drop = FALSE
    ]
  )

  x_test <- t(
    representation$new[
      genes,
      ,
      drop = FALSE
    ]
  )

  colnames(x_train) <- genes
  colnames(x_test) <- genes

  parameters <- architecture_lock_parameters(
    lock
  )

  seed_lgbm <- make_trace_subseed(
    GLOBAL_SEED,
    "stage04",
    "experiment4",
    disease_id,
    "final_LightGBM"
  )

  seed_tabpfn <- make_trace_subseed(
    GLOBAL_SEED,
    "stage04",
    "experiment4",
    disease_id,
    "final_TabPFN"
  )

  lgbm_model <- fit_python_lightgbm_object(
    x_train = x_train,
    y_train = split$development$metadata$outcome,
    parameters = parameters,
    seed = seed_lgbm,
    python = python
  )

  tabpfn_model <- fit_python_tabpfn_object(
    x_train = x_train,
    y_train = split$development$metadata$outcome,
    seed = seed_tabpfn,
    python = python
  )

  lgbm_probability <- predict_python_object(
    lgbm_model,
    x_test,
    python
  )

  tabpfn_probability <- predict_python_object(
    tabpfn_model,
    x_test,
    python
  )

  prediction_dt <- data.table::rbindlist(
    list(
      data.table(
        disease = disease_id,
        model = "LightGBM",
        sample_id = split$held_out$metadata$sample_id,
        cohort_group = as.character(
          split$held_out$metadata$cohort_group
        ),
        truth = as.character(
          split$held_out$metadata$outcome
        ),
        probability = lgbm_probability
      ),
      data.table(
        disease = disease_id,
        model = "TabPFN",
        sample_id = split$held_out$metadata$sample_id,
        cohort_group = as.character(
          split$held_out$metadata$cohort_group
        ),
        truth = as.character(
          split$held_out$metadata$outcome
        ),
        probability = tabpfn_probability
      )
    )
  )

  safe_write_csv(
    prediction_dt,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_architecture_frozen_heldout_predictions.csv"
      )
    )
  )

  final_metric_rows <- list()

  for (model_name in c(
    "LightGBM",
    "TabPFN"
  )) {
    probability <- if (
      model_name ==
        "LightGBM"
    ) {
      lgbm_probability
    } else {
      tabpfn_probability
    }

    bundle <- save_performance_bundle(
      output_dir = file.path(
        disease_dir,
        "frozen_heldout",
        model_name
      ),
      prefix = paste0(
        disease_id,
        "_",
        model_name,
        "_frozen_heldout"
      ),
      truth = split$held_out$metadata$outcome,
      probability = probability,
      threshold = 0.50,
      calibration_bins = 8L,
      bootstrap_iterations = 2000L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "stage04",
        "experiment4",
        disease_id,
        model_name,
        "final_curves"
      ),
      save_bootstrap = TRUE
    )

    current <- data.table::copy(
      bundle$metrics
    )

    current[
      ,
      `:=`(
        disease = disease_id,
        model = model_name
      )
    ]

    final_metric_rows[[
      length(final_metric_rows) + 1L
    ]] <- current
  }

  final_metrics <- data.table::rbindlist(
    final_metric_rows,
    use.names = TRUE,
    fill = TRUE
  )

  safe_write_csv(
    final_metrics,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_architecture_frozen_heldout_metrics.csv"
      )
    )
  )

  comparison <- paired_delong_test(
    truth = split$held_out$metadata$outcome,
    probability_a = lgbm_probability,
    probability_b = tabpfn_probability,
    label_a = "LightGBM",
    label_b = "TabPFN"
  )

  safe_write_csv(
    comparison,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_LightGBM_vs_TabPFN_paired_DeLong.csv"
      )
    )
  )

  model_dir <- file.path(
    disease_dir,
    "saved_models"
  )

  make_dir(
    model_dir
  )

  lgbm_pickle <- file.path(
    model_dir,
    paste0(
      disease_id,
      "_FINAL_LightGBM_model.pickle"
    )
  )

  tabpfn_pickle <- file.path(
    model_dir,
    paste0(
      disease_id,
      "_FINAL_TabPFN_model.pickle"
    )
  )

  reticulate::py_save_object(
    lgbm_model,
    lgbm_pickle
  )

  reticulate::py_save_object(
    tabpfn_model,
    tabpfn_pickle
  )

  x_train_path <- file.path(
    model_dir,
    paste0(
      disease_id,
      "_architecture_final_x_development.rds"
    )
  )

  x_test_path <- file.path(
    model_dir,
    paste0(
      disease_id,
      "_architecture_final_x_heldout.rds"
    )
  )

  safe_save_rds(
    x_train,
    x_train_path
  )

  safe_save_rds(
    x_test,
    x_test_path
  )

  state <- list(
    disease = disease_id,
    representation = "log2(CPM + 1)",
    feature_selection =
      "RF permutation importance",
    selected_genes = genes,
    top_n = 50L,
    LightGBM_config_id =
      as.character(
        lock$LightGBM_config_id[[1L]]
      ),
    LightGBM_parameters = parameters,
    positive_class = POSITIVE_CLASS,
    negative_class = NEGATIVE_CLASS,
    development_sample_ids =
      split$development$metadata$sample_id,
    heldout_sample_ids =
      split$held_out$metadata$sample_id,
    tabpfn_device = python$device,
    original_probabilities = prediction_dt
  )

  safe_save_rds(
    state,
    file.path(
      model_dir,
      paste0(
        disease_id,
        "_LightGBM_TabPFN_experiment_state.rds"
      )
    )
  )

  restored_lgbm <- reticulate::py_load_object(
    lgbm_pickle
  )

  restored_tabpfn <- reticulate::py_load_object(
    tabpfn_pickle
  )

  restored_lgbm_probability <- predict_python_object(
    restored_lgbm,
    x_test,
    python
  )

  restored_tabpfn_probability <- predict_python_object(
    restored_tabpfn,
    x_test,
    python
  )

  restoration_audit <- data.table(
    disease = disease_id,
    model = c(
      "LightGBM",
      "TabPFN"
    ),
    max_absolute_difference = c(
      max(
        abs(
          restored_lgbm_probability -
            lgbm_probability
        ),
        na.rm = TRUE
      ),
      max(
        abs(
          restored_tabpfn_probability -
            tabpfn_probability
        ),
        na.rm = TRUE
      )
    ),
    tolerance = 1e-5
  )

  restoration_audit[
    ,
    passed :=
      max_absolute_difference <=
        tolerance
  ]

  safe_write_csv(
    restoration_audit,
    file.path(
      disease_dir,
      paste0(
        disease_id,
        "_architecture_model_restoration_audit.csv"
      )
    )
  )

  if (!all(restoration_audit$passed)) {
    stop(
      disease_id,
      " restored architecture model probabilities exceeded the 1e-5 tolerance.",
      call. = FALSE
    )
  }

  # Interpretation begins only after restoration passes.
  molecular_dir <- file.path(
    disease_dir,
    "molecular_reliability"
  )

  make_dir(
    molecular_dir
  )

  lgbm_predict <- function(newdata) {
    predict_python_object(
      restored_lgbm,
      as.matrix(newdata),
      python
    )
  }

  tabpfn_predict <- function(newdata) {
    predict_python_object(
      restored_tabpfn,
      as.matrix(newdata),
      python
    )
  }

  lgbm_importance <- run_delta_auroc_permutation(
    model_label = paste0(
      disease_id,
      "_LightGBM"
    ),
    x = as.data.frame(
      x_test,
      check.names = FALSE
    ),
    truth = split$held_out$metadata$outcome,
    predict_probability = lgbm_predict,
    permutations = 100L,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "stage04",
      "experiment4",
      disease_id,
      "LightGBM",
      "importance"
    ),
    output_dir = file.path(
      molecular_dir,
      "LightGBM_importance"
    ),
    resume = TRUE
  )

  tabpfn_importance <- run_delta_auroc_permutation(
    model_label = paste0(
      disease_id,
      "_TabPFN"
    ),
    x = as.data.frame(
      x_test,
      check.names = FALSE
    ),
    truth = split$held_out$metadata$outcome,
    predict_probability = tabpfn_predict,
    permutations = 100L,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "stage04",
      "experiment4",
      disease_id,
      "TabPFN",
      "importance"
    ),
    output_dir = file.path(
      molecular_dir,
      "TabPFN_importance"
    ),
    resume = TRUE
  )

  concordance <- merge(
    lgbm_importance$summary[
      ,
      .(
        Gene,
        LightGBM_Delta_AUROC =
          Mean_Delta_AUROC,
        LightGBM_Rank = Rank
      )
    ],
    tabpfn_importance$summary[
      ,
      .(
        Gene,
        TabPFN_Delta_AUROC =
          Mean_Delta_AUROC,
        TabPFN_Rank = Rank
      )
    ],
    by = "Gene",
    all = TRUE
  )

  complete_concordance <- concordance[
    is.finite(
      LightGBM_Delta_AUROC
    ) &
      is.finite(
        TabPFN_Delta_AUROC
      )
  ]

  spearman <- if (
    nrow(complete_concordance) >= 3L
  ) {
    suppressWarnings(
      stats::cor(
        complete_concordance$LightGBM_Delta_AUROC,
        complete_concordance$TabPFN_Delta_AUROC,
        method = "spearman"
      )
    )
  } else {
    NA_real_
  }

  top_lgbm <- lgbm_importance$summary[
    order(-Mean_Delta_AUROC)
  ][
    seq_len(
      min(
        20L,
        .N
      )
    ),
    Gene
  ]

  top_tabpfn <- tabpfn_importance$summary[
    order(-Mean_Delta_AUROC)
  ][
    seq_len(
      min(
        20L,
        .N
      )
    ),
    Gene
  ]

  shared_top <- intersect(
    top_lgbm,
    top_tabpfn
  )

  union_top <- union(
    top_lgbm,
    top_tabpfn
  )

  jaccard <- if (length(union_top) > 0L) {
    length(shared_top) /
      length(union_top)
  } else {
    NA_real_
  }

  concordance_summary <- data.table(
    disease = disease_id,
    Spearman_all_selected_genes = spearman,
    Top_n = 20L,
    Shared_top_genes = length(shared_top),
    Union_top_genes = length(union_top),
    Top20_Jaccard = jaccard
  )

  safe_write_csv(
    concordance,
    file.path(
      molecular_dir,
      paste0(
        disease_id,
        "_cross_architecture_importance_alignment.csv"
      )
    )
  )

  safe_write_csv(
    concordance_summary,
    file.path(
      molecular_dir,
      paste0(
        disease_id,
        "_cross_architecture_concordance_summary.csv"
      )
    )
  )

  safe_write_csv(
    data.table(
      Gene = shared_top
    ),
    file.path(
      molecular_dir,
      paste0(
        disease_id,
        "_shared_top20_genes.csv"
      )
    )
  )

  scatter <- ggplot2::ggplot(
    complete_concordance,
    ggplot2::aes(
      x = LightGBM_Delta_AUROC,
      y = TabPFN_Delta_AUROC
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::labs(
      x = "LightGBM mean Delta-AUROC",
      y = "TabPFN mean Delta-AUROC",
      title = paste0(
        disease_id,
        " cross-architecture gene importance"
      )
    ) +
    ggplot2::theme_classic(
      base_size = 11
    )

  save_plot_dual(
    scatter,
    file.path(
      molecular_dir,
      paste0(
        disease_id,
        "_cross_architecture_importance_scatter"
      )
    ),
    width = 7,
    height = 6
  )

  if (
    requireNamespace(
      "DALEX",
      quietly = TRUE
    ) &&
      requireNamespace(
        "ingredients",
        quietly = TRUE
      )
  ) {
    lgbm_predict_function <- make_python_prediction_function(
      restored_lgbm,
      python
    )

    tabpfn_predict_function <- make_python_prediction_function(
      restored_tabpfn,
      python
    )

    lgbm_explainer <- make_dalex_explainer(
      model = restored_lgbm,
      x = as.data.frame(
        x_test,
        check.names = FALSE
      ),
      truth = split$held_out$metadata$outcome,
      predict_function = lgbm_predict_function,
      label = paste0(
        disease_id,
        "_LightGBM"
      )
    )

    tabpfn_explainer <- make_dalex_explainer(
      model = restored_tabpfn,
      x = as.data.frame(
        x_test,
        check.names = FALSE
      ),
      truth = split$held_out$metadata$outcome,
      predict_function = tabpfn_predict_function,
      label = paste0(
        disease_id,
        "_TabPFN"
      )
    )

    lgbm_ale_genes <- lgbm_importance$summary[
      order(-Mean_Delta_AUROC)
    ][
      seq_len(
        min(
          6L,
          .N
        )
      ),
      Gene
    ]

    tabpfn_ale_genes <- tabpfn_importance$summary[
      order(-Mean_Delta_AUROC)
    ][
      seq_len(
        min(
          6L,
          .N
        )
      ),
      Gene
    ]

    safe_write_csv(
      data.table(
        Model = c(
          rep(
            "LightGBM",
            length(lgbm_ale_genes)
          ),
          rep(
            "TabPFN",
            length(tabpfn_ale_genes)
          )
        ),
        Rank = c(
          seq_along(lgbm_ale_genes),
          seq_along(tabpfn_ale_genes)
        ),
        Gene = c(
          lgbm_ale_genes,
          tabpfn_ale_genes
        )
      ),
      file.path(
        molecular_dir,
        paste0(
          disease_id,
          "_architecture_ALE_selected_genes.csv"
        )
      )
    )

    run_ale_profiles(
      explainer = lgbm_explainer,
      genes = lgbm_ale_genes,
      output_dir = file.path(
        molecular_dir,
        "ALE"
      ),
      prefix = paste0(
        disease_id,
        "_LightGBM"
      )
    )

    run_ale_profiles(
      explainer = tabpfn_explainer,
      genes = tabpfn_ale_genes,
      output_dir = file.path(
        molecular_dir,
        "ALE"
      ),
      prefix = paste0(
        disease_id,
        "_TabPFN"
      )
    )
  } else {
    warning(
      disease_id,
      ": DALEX/ingredients unavailable; architecture ALE skipped."
    )
  }

  list(
    disease = disease_id,
    loco = loco,
    final_metrics = final_metrics,
    final_predictions = prediction_dt,
    paired_delong = comparison,
    restoration_audit = restoration_audit,
    concordance = concordance_summary
  )
}

run_experiment_4_testing <- function() {
  message_section(
    "Stage 04 | Experiment 4 | model architecture"
  )

  output_root <- file.path(
    STAGE4_ROOT,
    "04_experiment_4_model_architecture"
  )

  make_dir(
    output_root
  )

  require_packages(
    "reticulate",
    bioc = FALSE
  )

  python <- stage3$initialise_architecture_python()

  results <- list()

  for (disease_id in REQUESTED_DISEASES) {
    results[[disease_id]] <- run_experiment_4_final_disease(
      disease_id = disease_id,
      python = python,
      output_root = output_root
    )
  }

  metrics <- data.table::rbindlist(
    lapply(
      results,
      function(x) x$final_metrics
    ),
    use.names = TRUE,
    fill = TRUE
  )

  loco <- data.table::rbindlist(
    lapply(
      results,
      function(x) {
        current <- data.table::copy(
          x$loco$summary
        )

        current[
          ,
          disease := x$disease
        ]

        current
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  comparisons <- data.table::rbindlist(
    lapply(
      results,
      function(x) {
        current <- data.table::copy(
          x$paired_delong
        )

        current[
          ,
          disease := x$disease
        ]

        current
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  concordance <- data.table::rbindlist(
    lapply(
      results,
      function(x) x$concordance
    ),
    use.names = TRUE,
    fill = TRUE
  )

  safe_write_csv(
    metrics,
    file.path(
      output_root,
      "Experiment4_all_disease_frozen_heldout_metrics.csv"
    )
  )

  safe_write_csv(
    loco,
    file.path(
      output_root,
      "Experiment4_all_disease_LOCO_summary.csv"
    )
  )

  safe_write_csv(
    comparisons,
    file.path(
      output_root,
      "Experiment4_all_disease_LightGBM_vs_TabPFN_DeLong.csv"
    )
  )

  safe_write_csv(
    concordance,
    file.path(
      output_root,
      "Experiment4_all_disease_cross_architecture_concordance.csv"
    )
  )

  invisible(
    list(
      metrics = metrics,
      loco = loco,
      comparisons = comparisons,
      concordance = concordance
    )
  )
}



# ==============================================================================
# 11. STAGE-04 METHOD MANIFEST
# ==============================================================================

method_manifest <- data.table(
  experiment = c(
    "Experiment 1",
    "Experiment 2",
    "Experiment 3",
    "Experiment 4"
  ),
  primary_testing_design = c(
    "Frozen balanced 16-sample held-out evaluation",
    "LODO plus identical frozen TB comparison",
    "Selected development LOCO plus frozen disease-specific held-out evaluation",
    "Matched architecture LOCO plus identical frozen disease-specific held-out evaluation"
  ),
  statistical_reliability = c(
    "Exact permutation; paired DeLong; 5000 gap bootstrap",
    "Paired DeLong; 2000 ROC bootstrap",
    "DeLong AUROC CI; 2000 ROC/PR bootstrap; calibration",
    "Paired DeLong; 2000 ROC/PR bootstrap; calibration; restoration audit"
  ),
  molecular_reliability = c(
    "Final top-50 gene vectors retained",
    "50-permutation Delta-AUROC plus Reactome",
    "100-permutation Delta-AUROC plus ALE",
    "100-permutation Delta-AUROC; Spearman; top20 Jaccard; ALE"
  ),
  held_out_used_for_model_selection = FALSE
)

safe_write_csv(
  method_manifest,
  file.path(
    AUDIT_ROOT,
    "TRACE_Omics_stage04_method_manifest.csv"
  )
)


# ==============================================================================
# 12. RUN REQUESTED EXPERIMENTS
# ==============================================================================

stage04_results <- list()

if (run_requested(
  "experiment1"
)) {
  stage04_results$Experiment1 <-
    run_experiment_1_testing()
}

if (run_requested(
  "experiment2"
)) {
  stage04_results$Experiment2 <-
    run_experiment_2_testing()
}

if (run_requested(
  "experiment3"
)) {
  stage04_results$Experiment3 <-
    run_experiment_3_testing()
}

if (run_requested(
  "experiment4"
)) {
  stage04_results$Experiment4 <-
    run_experiment_4_testing()
}


# ==============================================================================
# 13. FINAL REPRODUCIBILITY RECORDS
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

run_manifest <- data.table(
  stage = "04_model_testing_reliability",
  created_at = as.character(
    Sys.time()
  ),
  configuration = CONFIG_PATH,
  modes = paste(
    RUN_MODES,
    collapse = ","
  ),
  requested_diseases = paste(
    REQUESTED_DISEASES,
    collapse = ","
  ),
  global_seed = GLOBAL_SEED,
  strict_audits = STRICT_AUDITS,
  frozen_evaluation_used_for_selection = FALSE,
  experiment4_restore_tolerance = 1e-5
)

safe_write_csv(
  run_manifest,
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage04_run_manifest.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage04_sessionInfo.txt"
  )
)

safe_save_rds(
  stage04_results,
  file.path(
    STAGE4_ROOT,
    "TRACE_Omics_stage04_master_results.rds"
  )
)

writeLines(
  c(
    "TRACE-Omics Stage 04 - Model testing and reliability assessment",
    "",
    "This stage consumes Stage-03 development-only locks.",
    "Frozen held-out populations are used for final scoring and post-hoc reliability analysis.",
    "Frozen held-out outcomes are never used to revise representation, features, class handling, learner family, hyperparameters or operating threshold.",
    "",
    "Experiment 1:",
    "  equal-size frozen testing + exact permutation + paired DeLong + 5,000 gap bootstrap",
    "",
    "Experiment 2:",
    "  LODO + frozen Model A versus Model B TB comparison + permutation importance + Reactome",
    "",
    "Experiment 3:",
    "  selected LOCO transportability + frozen disease-specific Random Forest + permutation importance + ALE",
    "",
    "Experiment 4:",
    "  matched LightGBM/TabPFN LOCO + frozen final comparison + model restoration + cross-architecture molecular reliability",
    "",
    "For Experiment 4, execute one disease per fresh R session when disease-specific Python/Conda environments are used."
  ),
  con = file.path(
    STAGE4_ROOT,
    "README_STAGE04.txt"
  )
)


# ==============================================================================
# 14. FINAL CONSOLE OUTPUT
# ==============================================================================

cat(
  "\n\n============================================================\n",
  "TRACE-OMICS STAGE 04 COMPLETE\n",
  "MODEL TESTING AND RELIABILITY ASSESSMENT\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nOutput root:\n",
  STAGE4_ROOT,
  "\n",
  sep = ""
)

cat(
  "\nFrozen evaluation used for model selection: FALSE\n",
  sep = ""
)

cat(
  "\nCompleted modes: ",
  paste(
    RUN_MODES,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "STAGE 04 REPRODUCIBILITY RECORDS WRITTEN\n",
  "============================================================\n",
  sep = ""
)

