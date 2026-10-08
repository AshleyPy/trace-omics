# RELEASE NOTE -----------------------------------------------------------------
# This source is a recovered comprehensive TRACE-Omics implementation from the
# dissertation code-development history. The final statistical design is defined
# in config/final_dissertation.yml and docs/METHODS_TO_CODE.md. Earlier branches
# retained inside this file are provenance-bearing and may differ from the final
# reported sensitivity estimand. See docs/REPRODUCTION_STATUS.md before rerunning.
# -----------------------------------------------------------------------------
# ==============================================================================
# TRACE-Omics
# Stage 02: Feature-set characterisation
# ==============================================================================
#
# File:
#   scripts/02_feature_set_characterisation/02_feature_set_characterisation.R
#
# Framework:
#   Transcriptomic Robustness Assessment and Cross-cohort Evaluation for
#   Omics learning (TRACE-Omics)
#
# Purpose:
#   Reconstruct the complete feature-set characterisation stage used by the
#   dissertation after Stage 01 representation specification and before model
#   development/optimisation.
#
#   The script:
#     1. loads frozen development/source representation matrices from Stage 01;
#     2. applies the eight full-development feature-selection procedures;
#     3. stores complete feature-score tables and deterministic ranked gene lists;
#     4. freezes top-N gene sets using fitting data only;
#     5. quantifies feature-set concordance with Jaccard similarity;
#     6. quantifies redundancy/conditioning using VIF and effective rank;
#     7. bootstraps conditioning metrics while keeping each selected gene set fixed;
#     8. prepares experiment-specific feature artefacts for:
#          Experiment 1 - sample-support perturbation;
#          Experiment 2 - disease-composition / LODO perturbation;
#          Experiment 3 - cohort-composition / LOCO perturbation;
#          Experiment 4 - matched LightGBM-TabPFN architecture perturbation;
#     9. saves fold/source-specific selected-gene vectors, score tables, selected
#        matrices, manifests, diagnostics and session information.
#
# Leakage control:
#   - Feature ranking is learned only from the fitting partition represented by
#     each Stage 01 development/source matrix.
#   - Evaluation/held-out/LODO/LOCO observations are only subset to the frozen
#     gene vector after selection.
#   - Outer-fold architecture gene selection is repeated inside each outer
#     training fold.
#   - Nested-CV inner-fold feature selection remains a Stage 03 responsibility;
#     Stage 02 does not pre-select one global set for use inside inner tuning.
#
# Full-development selectors:
#   ANOVA
#   AUC
#   RF permutation importance
#   Boruta
#   CMIM
#   mRMR
#   Information gain
#   Variance
#
# Canonical full-development top-N:
#   50 genes
#
# Experiment-specific feature rules:
#   Experiment 1:
#     Information gain, RF importance, Boruta, mRMR, CMIM; top 50.
#
#   Experiment 2:
#     Model A/B: mRMR and CMIM; top 50.
#     Locked LODO: mRMR; top 50.
#
#   Experiment 3:
#     RF permutation importance; candidate top 25 / 50 / 100.
#
#   Experiment 4:
#     RF permutation importance; top 50; matched across LightGBM and TabPFN.
#
# Example execution:
#
#   Rscript scripts/02_feature_set_characterisation/02_feature_set_characterisation.R
#
#   Rscript scripts/02_feature_set_characterisation/02_feature_set_characterisation.R \
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

    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]

    key <- pieces[[1]]
    value <- if (length(pieces) >= 2L) paste(pieces[-1], collapse = "=") else "true"

    output[[key]] <- value
  }

  output
}

parse_boolean <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }

  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

args <- parse_named_arguments(commandArgs(trailingOnly = TRUE))

CONFIG_PATH <- if (!is.null(args$config)) args$config else "config/example_config.yml"
MODE_TEXT <- if (!is.null(args$mode)) args$mode else "all"
DISEASE_TEXT <- if (!is.null(args$diseases)) args$diseases else "TB,RA,UC,CD,T1D"
STRICT_AUDITS <- parse_boolean(args$strict, default = FALSE)

RUN_MODES <- unique(trimws(strsplit(MODE_TEXT, ",", fixed = TRUE)[[1]]))
REQUESTED_DISEASES <- unique(trimws(strsplit(DISEASE_TEXT, ",", fixed = TRUE)[[1]]))

VALID_MODES <- c(
  "all",
  "full",
  "experiment1",
  "experiment2",
  "experiment3",
  "experiment4"
)

if (!all(RUN_MODES %in% VALID_MODES)) {
  stop(
    "Unknown --mode value(s): ",
    paste(setdiff(RUN_MODES, VALID_MODES), collapse = ", "),
    call. = FALSE
  )
}

run_requested <- function(mode) {
  "all" %in% RUN_MODES || mode %in% RUN_MODES
}


# ==============================================================================
# 1. LOAD TRACE-OMICS UTILITIES, PACKAGES AND CONFIGURATION
# ==============================================================================

UTILS_PATH <- file.path("R", "utils_io.R")

if (!file.exists(UTILS_PATH)) {
  stop(
    "TRACE-Omics utility file not found: ",
    UTILS_PATH,
    "\nRun this script from the repository root.",
    call. = FALSE
  )
}

source(UTILS_PATH)

require_packages(
  c(
    "yaml",
    "data.table",
    "ggplot2",
    "matrixStats",
    "mlr3",
    "mlr3filters",
    "mlr3learners",
    "ranger",
    "praznik",
    "FSelectorRcpp",
    "Boruta",
    "pROC"
  ),
  bioc = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

config <- read_trace_config(CONFIG_PATH)

GLOBAL_SEED <- as.integer(get_config_value(config, "project.seed", 7832))
POSITIVE_CLASS <- as.character(get_config_value(config, "project.positive_class", "case"))
NEGATIVE_CLASS <- as.character(get_config_value(config, "project.negative_class", "control"))

DEFAULT_TOP_N <- as.integer(
  get_config_value(
    config,
    "feature_selection.default_top_n",
    get_config_value(config, "analysis.top_n_features", 50)
  )
)

CONDITIONING_BOOTSTRAPS <- as.integer(
  get_config_value(
    config,
    "feature_set_characterisation.conditioning.bootstrap_iterations",
    200
  )
)

CONDITIONING_CONFIDENCE <- as.numeric(
  get_config_value(
    config,
    "feature_set_characterisation.conditioning.confidence_level",
    0.95
  )
)

MINIMUM_FINITE_BOOTSTRAPS <- as.integer(
  get_config_value(
    config,
    "feature_set_characterisation.conditioning.minimum_finite_bootstraps",
    20
  )
)

RF_IMPORTANCE_TREES_FULL <- 300L

RF_IMPORTANCE_TREES_DISEASE_SPECIFIC <- as.integer(
  get_config_value(
    config,
    "cohort_composition_experiment.final_learner.rf_importance_trees",
    1000
  )
)

RF_IMPORTANCE_TREES_ARCHITECTURE <- as.integer(
  get_config_value(
    config,
    "model_architecture_experiment.feature_selection_trees",
    RF_IMPORTANCE_TREES_DISEASE_SPECIFIC
  )
)

set_trace_seed(GLOBAL_SEED)


# ==============================================================================
# 2. FEATURE-SELECTION REGISTRY
# ==============================================================================

FULL_SELECTOR_IDS <- c(
  "anova",
  "auc",
  "importance",
  "boruta",
  "cmim",
  "mrmr",
  "information_gain",
  "variance"
)

FULL_SELECTOR_LABELS <- c(
  anova = "ANOVA",
  auc = "AUC",
  importance = "RF importance",
  boruta = "Boruta",
  cmim = "CMIM",
  mrmr = "mRMR",
  information_gain = "Information gain",
  variance = "Variance"
)

EXPERIMENT1_SELECTOR_IDS <- c(
  "information_gain",
  "importance",
  "boruta",
  "mrmr",
  "cmim"
)

EXPERIMENT2_MODEL_SELECTOR_IDS <- c(
  "mrmr",
  "cmim"
)

EXPERIMENT2_LOCKED_SELECTOR_ID <- "mrmr"

EXPERIMENT3_SELECTOR_ID <- "importance"
EXPERIMENT3_FEATURE_COUNTS <- c(25L, 50L, 100L)

EXPERIMENT4_SELECTOR_ID <- "importance"
EXPERIMENT4_TOP_N <- as.integer(
  get_config_value(
    config,
    "model_architecture_experiment.top_n_features",
    50
  )
)

REPRESENTATION_LABELS <- c(
  RAW_counts = "Raw counts",
  VST_only = "VST only",
  VST_ComBat = "VST + ComBat",
  VST_ComBatRef = "VST + ComBat-Ref",
  VST_fSVA = "VST + fSVA",
  ComBatSeq_log1p = "ComBat-Seq + log1p",
  log2_CPM_plus_1 = "log2(CPM + 1)"
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
    paste(setdiff(REQUESTED_DISEASES, ALL_DISEASES), collapse = ", "),
    call. = FALSE
  )
}


# ==============================================================================
# 3. STAGE INPUT/OUTPUT ROOTS
# ==============================================================================

stage1_example <- get_config_value(
  config,
  "paths.stage_results.representation_specification",
  "results/01_representation_specification/TB"
)

STAGE1_ROOT <- dirname(stage1_example)

stage2_example <- get_config_value(
  config,
  "paths.stage_results.feature_set_characterisation",
  "results/02_feature_set_characterisation/TB"
)

STAGE2_ROOT <- dirname(stage2_example)

AUDIT_ROOT <- file.path(STAGE2_ROOT, "00_audits")
SESSION_ROOT <- file.path(STAGE2_ROOT, "99_session_info")

make_dir(STAGE2_ROOT)
make_dir(AUDIT_ROOT)
make_dir(SESSION_ROOT)

message_section("TRACE-Omics Stage 02: Feature-set characterisation")
message("Configuration:          ", CONFIG_PATH)
message("Modes:                  ", paste(RUN_MODES, collapse = ", "))
message("Diseases:               ", paste(REQUESTED_DISEASES, collapse = ", "))
message("Global seed:            ", GLOBAL_SEED)
message("Top-N:                  ", DEFAULT_TOP_N)
message("Conditioning bootstraps:", CONDITIONING_BOOTSTRAPS)
message("Strict audits:          ", STRICT_AUDITS)


# ==============================================================================
# 4. GENERAL FILE/DATA HELPERS
# ==============================================================================

first_matching_file <- function(directory,
                                pattern,
                                required = TRUE,
                                label = "file") {
  if (!dir.exists(directory)) {
    if (isTRUE(required)) {
      stop(
        "Directory not found while locating ",
        label,
        ": ",
        directory,
        call. = FALSE
      )
    }

    return(NA_character_)
  }

  matches <- sort(
    list.files(
      directory,
      pattern = pattern,
      full.names = TRUE,
      recursive = FALSE
    )
  )

  if (length(matches) == 0L) {
    if (isTRUE(required)) {
      stop(
        "No ",
        label,
        " matching ",
        pattern,
        " found in ",
        directory,
        call. = FALSE
      )
    }

    return(NA_character_)
  }

  matches[[1]]
}

read_matrix_rds <- function(path,
                            object_name = "matrix") {
  mat <- safe_read_rds(path, object_name)
  mat <- coerce_numeric_matrix(mat, object_name)

  if (is.null(rownames(mat)) || anyDuplicated(rownames(mat)) > 0L) {
    stop(
      object_name,
      " requires unique gene row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(mat)) || anyDuplicated(colnames(mat)) > 0L) {
    stop(
      object_name,
      " requires unique sample column names.",
      call. = FALSE
    )
  }

  mat
}

read_metadata_csv <- function(path,
                              object_name = "metadata") {
  metadata <- as.data.frame(
    safe_read_csv(path),
    check.names = FALSE
  )

  check_required_columns(
    metadata,
    c("sample_id", "outcome"),
    object_name
  )

  metadata$sample_id <- trimws(as.character(metadata$sample_id))

  if ("cohort_group" %in% colnames(metadata)) {
    metadata$cohort_group <- factor(as.character(metadata$cohort_group))
  } else if ("cohort" %in% colnames(metadata)) {
    metadata$cohort_group <- factor(as.character(metadata$cohort))
  }

  metadata$outcome <- factor(
    as.character(metadata$outcome),
    levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)
  )

  if (anyDuplicated(metadata$sample_id) > 0L) {
    stop(
      object_name,
      " contains duplicated sample IDs.",
      call. = FALSE
    )
  }

  assert_binary_outcome(
    metadata,
    outcome_col = "outcome",
    negative_class = NEGATIVE_CLASS,
    positive_class = POSITIVE_CLASS
  )

  metadata
}

align_matrix_metadata <- function(mat,
                                  metadata,
                                  object_name = "feature-selection input") {
  shared <- intersect(
    colnames(mat),
    metadata$sample_id
  )

  if (length(shared) < 2L) {
    stop(
      object_name,
      " has fewer than two aligned observations.",
      call. = FALSE
    )
  }

  mat <- mat[
    ,
    shared,
    drop = FALSE
  ]

  metadata <- metadata[
    match(shared, metadata$sample_id),
    ,
    drop = FALSE
  ]

  if (!identical(colnames(mat), metadata$sample_id)) {
    stop(
      object_name,
      " sample alignment failed.",
      call. = FALSE
    )
  }

  list(
    matrix = mat,
    metadata = metadata
  )
}

normalise_representation_label <- function(x) {
  if (x %in% names(REPRESENTATION_LABELS)) {
    return(unname(REPRESENTATION_LABELS[[x]]))
  }

  x
}

normalise_selector_label <- function(x) {
  if (x %in% names(FULL_SELECTOR_LABELS)) {
    return(unname(FULL_SELECTOR_LABELS[[x]]))
  }

  x
}

deterministic_rank_table <- function(score_table) {
  score_table <- as.data.table(score_table)

  score_table <- score_table[
    !is.na(gene) &
      nzchar(gene) &
      is.finite(score)
  ]

  setorder(
    score_table,
    -score,
    gene
  )

  score_table[
    ,
    rank := seq_len(.N)
  ]

  score_table
}

save_gene_vector <- function(genes,
                             path) {
  make_dir(dirname(path))

  genes <- unique(
    as.character(genes)
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  safe_save_rds(
    genes,
    path
  )

  invisible(genes)
}

save_selected_matrix_pair <- function(train_matrix,
                                      evaluation_matrix,
                                      genes,
                                      output_dir,
                                      prefix) {
  make_dir(output_dir)

  genes <- unique(as.character(genes))

  missing_train <- setdiff(
    genes,
    rownames(train_matrix)
  )

  if (length(missing_train) > 0L) {
    stop(
      "Selected genes missing from fitting matrix: ",
      paste(missing_train, collapse = ", "),
      call. = FALSE
    )
  }

  selected_train <- train_matrix[
    genes,
    ,
    drop = FALSE
  ]

  safe_save_rds(
    selected_train,
    file.path(
      output_dir,
      paste0(prefix, "_development_selected_matrix.rds")
    )
  )

  if (!is.null(evaluation_matrix)) {
    missing_eval <- setdiff(
      genes,
      rownames(evaluation_matrix)
    )

    if (length(missing_eval) > 0L) {
      stop(
        "Selected genes missing from evaluation matrix: ",
        paste(missing_eval, collapse = ", "),
        call. = FALSE
      )
    }

    selected_evaluation <- evaluation_matrix[
      genes,
      ,
      drop = FALSE
    ]

    assert_feature_alignment(
      selected_train,
      selected_evaluation,
      reference_label = "selected fitting matrix",
      target_label = "selected evaluation matrix"
    )

    safe_save_rds(
      selected_evaluation,
      file.path(
        output_dir,
        paste0(prefix, "_evaluation_selected_matrix.rds")
      )
    )
  }

  invisible(selected_train)
}


# ==============================================================================
# 5. MLR3 FEATURE-TASK CONSTRUCTION
# ==============================================================================

make_safe_gene_names <- function(genes) {
  make.unique(
    make.names(
      genes,
      unique = TRUE
    )
  )
}

make_feature_task <- function(expression_matrix,
                              metadata,
                              task_id) {
  expression_matrix <- coerce_numeric_matrix(
    expression_matrix,
    paste0(task_id, " expression matrix")
  )

  aligned <- align_matrix_metadata(
    expression_matrix,
    metadata,
    paste0(task_id, " feature-selection input")
  )

  expression_matrix <- aligned$matrix
  metadata <- aligned$metadata

  original_genes <- rownames(expression_matrix)
  safe_genes <- make_safe_gene_names(original_genes)

  x <- as.data.table(
    t(expression_matrix)
  )

  setnames(
    x,
    safe_genes
  )

  x[
    ,
    outcome := factor(
      as.character(metadata$outcome),
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


# ==============================================================================
# 6. FEATURE FILTER OBJECTS
# ==============================================================================

make_filter_object <- function(filter_id,
                               seed,
                               rf_trees = RF_IMPORTANCE_TREES_FULL) {
  if (filter_id == "importance") {
    importance_learner <- mlr3::lrn(
      "classif.ranger",
      predict_type = "prob",
      importance = "permutation",
      num.trees = as.integer(rf_trees),
      num.threads = 1L,
      seed = as.integer(seed)
    )

    return(
      mlr3filters::flt(
        "importance",
        learner = importance_learner
      )
    )
  }

  mlr3filters::flt(
    filter_id
  )
}


# ==============================================================================
# 7. TRAINING-ONLY FEATURE SELECTION
# ==============================================================================

select_ranked_features <- function(train_matrix,
                                   train_metadata,
                                   filter_id,
                                   top_n = DEFAULT_TOP_N,
                                   seed = GLOBAL_SEED,
                                   rf_trees = RF_IMPORTANCE_TREES_FULL) {
  if (!filter_id %in% FULL_SELECTOR_IDS) {
    stop(
      "Unknown feature-selection method: ",
      filter_id,
      call. = FALSE
    )
  }

  train_matrix <- coerce_numeric_matrix(
    train_matrix,
    paste0(filter_id, " fitting matrix")
  )

  training_variance <- matrixStats::rowVars(
    train_matrix
  )

  keep <- (
    is.finite(training_variance) &
      training_variance > sqrt(.Machine$double.eps)
  )

  train_matrix <- train_matrix[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(train_matrix) < top_n) {
    stop(
      filter_id,
      ": fewer than ",
      top_n,
      " variable genes were available.",
      call. = FALSE
    )
  }

  task_info <- make_feature_task(
    expression_matrix = train_matrix,
    metadata = train_metadata,
    task_id = paste0(
      "TRACE_FS_",
      safe_file_name(filter_id)
    )
  )

  filter_object <- make_filter_object(
    filter_id = filter_id,
    seed = seed,
    rf_trees = rf_trees
  )

  set.seed(seed)

  filter_object$calculate(
    task_info$task
  )

  scores <- filter_object$scores

  score_table <- data.table(
    safe_gene = names(scores),
    score = as.numeric(scores)
  )

  score_table <- score_table[
    !is.na(safe_gene) &
      nzchar(safe_gene) &
      is.finite(score)
  ]

  score_table[
    ,
    gene := task_info$map$original_gene[
      match(
        safe_gene,
        task_info$map$safe_gene
      )
    ]
  ]

  score_table <- score_table[
    !is.na(gene) &
      nzchar(gene)
  ]

  score_table <- deterministic_rank_table(
    score_table[
      ,
      .(
        gene,
        score
      )
    ]
  )

  if (nrow(score_table) < top_n) {
    stop(
      filter_id,
      ": fewer than ",
      top_n,
      " finite feature scores were returned.",
      call. = FALSE
    )
  }

  selected_genes <- score_table[
    seq_len(top_n),
    gene
  ]

  if (anyDuplicated(selected_genes) > 0L) {
    stop(
      filter_id,
      ": duplicated genes occurred in the selected feature vector.",
      call. = FALSE
    )
  }

  list(
    filter_id = filter_id,
    genes = selected_genes,
    scores = score_table,
    variable_genes = rownames(train_matrix)
  )
}

run_selector_collection <- function(train_matrix,
                                    train_metadata,
                                    selector_ids,
                                    top_n,
                                    seed,
                                    rf_trees = RF_IMPORTANCE_TREES_FULL,
                                    output_dir = NULL,
                                    prefix = NULL) {
  outputs <- list()

  for (selector_index in seq_along(selector_ids)) {
    selector_id <- selector_ids[[selector_index]]

    current_seed <- make_trace_subseed(
      seed,
      selector_id,
      selector_index
    )

    message(
      "  [feature selection] ",
      normalise_selector_label(selector_id),
      " | top ",
      top_n
    )

    result <- select_ranked_features(
      train_matrix = train_matrix,
      train_metadata = train_metadata,
      filter_id = selector_id,
      top_n = top_n,
      seed = current_seed,
      rf_trees = rf_trees
    )

    outputs[[selector_id]] <- result

    if (!is.null(output_dir)) {
      selector_dir <- file.path(
        output_dir,
        safe_file_name(selector_id)
      )

      make_dir(selector_dir)

      safe_write_csv(
        result$scores,
        file.path(
          selector_dir,
          paste0(
            prefix,
            "_",
            selector_id,
            "_full_ranking.csv"
          )
        )
      )

      save_gene_vector(
        result$genes,
        file.path(
          selector_dir,
          paste0(
            prefix,
            "_",
            selector_id,
            "_top",
            top_n,
            "_genes.rds"
          )
        )
      )

      safe_write_csv(
        data.table(
          rank = seq_along(result$genes),
          gene = result$genes
        ),
        file.path(
          selector_dir,
          paste0(
            prefix,
            "_",
            selector_id,
            "_top",
            top_n,
            "_genes.csv"
          )
        )
      )
    }
  }

  outputs
}


# ==============================================================================
# 8. JACCARD FEATURE-SET CONCORDANCE
# ==============================================================================

jaccard_index <- function(a,
                          b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))

  union_n <- length(
    union(a, b)
  )

  if (union_n == 0L) {
    return(NA_real_)
  }

  length(intersect(a, b)) /
    union_n
}

jaccard_matrix <- function(feature_sets) {
  method_names <- names(feature_sets)

  output <- matrix(
    NA_real_,
    nrow = length(method_names),
    ncol = length(method_names),
    dimnames = list(
      method_names,
      method_names
    )
  )

  for (i in seq_along(method_names)) {
    for (j in seq_along(method_names)) {
      output[i, j] <- jaccard_index(
        feature_sets[[method_names[[i]]]],
        feature_sets[[method_names[[j]]]]
      )
    }
  }

  diag(output) <- 1

  output
}

tidy_jaccard_matrix <- function(matrix_object,
                                disease,
                                representation,
                                top_n,
                                context) {
  dt <- as.data.table(
    as.table(matrix_object)
  )

  setnames(
    dt,
    c(
      "method_1",
      "method_2",
      "jaccard"
    )
  )

  dt[
    ,
    `:=`(
      disease = disease,
      representation = representation,
      top_n = top_n,
      context = context
    )
  ]

  dt
}

summarise_jaccard <- function(jaccard_dt) {
  off_diagonal <- jaccard_dt[
    method_1 != method_2
  ]

  if (nrow(off_diagonal) == 0L) {
    return(data.table())
  }

  representation_summary <- off_diagonal[
    ,
    .(
      mean_jaccard = mean(jaccard, na.rm = TRUE),
      median_jaccard = median(jaccard, na.rm = TRUE),
      min_jaccard = min(jaccard, na.rm = TRUE),
      max_jaccard = max(jaccard, na.rm = TRUE)
    ),
    by = .(
      disease,
      representation,
      top_n,
      context
    )
  ]

  method_summary <- off_diagonal[
    ,
    .(
      mean_method_jaccard = mean(jaccard, na.rm = TRUE),
      median_method_jaccard = median(jaccard, na.rm = TRUE)
    ),
    by = .(
      disease,
      representation,
      method = method_1,
      top_n,
      context
    )
  ]

  list(
    representation = representation_summary,
    method = method_summary
  )
}

plot_jaccard_heatmap <- function(jaccard_dt,
                                 output_path,
                                 title) {
  plot_dt <- copy(jaccard_dt)

  plot_dt[
    ,
    label := ifelse(
      is.finite(jaccard),
      sprintf("%.2f", jaccard),
      ""
    )
  ]

  p <- ggplot(
    plot_dt,
    aes(
      x = method_2,
      y = method_1,
      fill = jaccard
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.25
    ) +
    geom_text(
      aes(label = label),
      size = 2.6
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 1),
      na.value = "grey95",
      name = "Jaccard"
    ) +
    coord_equal() +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    output_path,
    p,
    width = 170,
    height = 150,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  invisible(p)
}


# ==============================================================================
# 9. FEATURE-SPACE CONDITIONING
# ==============================================================================

prepare_conditioning_matrix <- function(expression_matrix,
                                        genes) {
  genes <- intersect(
    unique(as.character(genes)),
    rownames(expression_matrix)
  )

  if (length(genes) < 2L) {
    return(NULL)
  }

  x <- t(
    expression_matrix[
      genes,
      ,
      drop = FALSE
    ]
  )

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  finite_columns <- colSums(!is.finite(x)) == 0L
  x <- x[, finite_columns, drop = FALSE]

  if (ncol(x) < 2L) {
    return(NULL)
  }

  column_sd <- apply(
    x,
    2L,
    stats::sd
  )

  keep <- (
    is.finite(column_sd) &
      column_sd > sqrt(.Machine$double.eps)
  )

  x <- x[, keep, drop = FALSE]

  if (ncol(x) < 2L) {
    return(NULL)
  }

  scale(x)
}

effective_rank <- function(x) {
  if (is.null(x) || ncol(x) < 2L) {
    return(NA_real_)
  }

  singular_values <- svd(
    x,
    nu = 0,
    nv = 0
  )$d

  singular_values <- singular_values[
    is.finite(singular_values) &
      singular_values > 0
  ]

  if (length(singular_values) == 0L) {
    return(NA_real_)
  }

  probabilities <- singular_values^2 /
    sum(singular_values^2)

  entropy <- -sum(
    probabilities *
      log(probabilities)
  )

  exp(entropy)
}

vif_vector <- function(x) {
  if (is.null(x) || ncol(x) < 2L) {
    return(numeric(0))
  }

  x <- as.matrix(x)
  p <- ncol(x)

  output <- rep(
    NA_real_,
    p
  )

  names(output) <- colnames(x)

  for (j in seq_len(p)) {
    y <- x[, j]

    predictors <- x[
      ,
      -j,
      drop = FALSE
    ]

    if (ncol(predictors) == 0L) {
      output[[j]] <- 1
      next
    }

    model <- tryCatch(
      stats::lm.fit(
        x = cbind(
          intercept = 1,
          predictors
        ),
        y = y
      ),
      error = function(e) NULL
    )

    if (is.null(model)) {
      next
    }

    fitted <- model$fitted.values

    sst <- sum(
      (y - mean(y))^2
    )

    sse <- sum(
      (y - fitted)^2
    )

    if (
      !is.finite(sst) ||
        sst <= sqrt(.Machine$double.eps)
    ) {
      next
    }

    r_squared <- 1 - sse / sst
    r_squared <- min(
      max(r_squared, 0),
      1
    )

    output[[j]] <- if (r_squared >= 1 - 1e-12) {
      Inf
    } else {
      1 / (1 - r_squared)
    }
  }

  output
}

conditioning_metrics_once <- function(expression_matrix,
                                      genes) {
  x <- prepare_conditioning_matrix(
    expression_matrix,
    genes
  )

  if (is.null(x)) {
    return(
      data.table(
        n_selected = length(unique(genes)),
        n_used = NA_integer_,
        mean_vif = NA_real_,
        max_vif = NA_real_,
        effective_rank = NA_real_,
        normalised_effective_rank = NA_real_,
        mean_abs_correlation = NA_real_
      )
    )
  }

  vifs <- vif_vector(x)

  finite_vifs <- vifs[
    is.finite(vifs)
  ]

  correlation_matrix <- suppressWarnings(
    stats::cor(x)
  )

  mean_abs_correlation <- if (ncol(x) >= 2L) {
    mean(
      abs(
        correlation_matrix[
          lower.tri(
            correlation_matrix
          )
        ]
      ),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  erank <- effective_rank(x)

  data.table(
    n_selected = length(unique(genes)),
    n_used = ncol(x),
    mean_vif = if (length(finite_vifs) > 0L) mean(finite_vifs) else Inf,
    max_vif = if (length(vifs) > 0L) max(vifs, na.rm = TRUE) else NA_real_,
    effective_rank = erank,
    normalised_effective_rank = erank / ncol(x),
    mean_abs_correlation = mean_abs_correlation
  )
}

bootstrap_conditioning_metrics <- function(expression_matrix,
                                           genes,
                                           n_boot = CONDITIONING_BOOTSTRAPS,
                                           seed = GLOBAL_SEED,
                                           confidence = CONDITIONING_CONFIDENCE) {
  point <- conditioning_metrics_once(
    expression_matrix,
    genes
  )

  boot_rows <- vector(
    "list",
    n_boot
  )

  for (bootstrap_index in seq_len(n_boot)) {
    set_trace_subseed(
      seed,
      "conditioning_bootstrap",
      bootstrap_index
    )

    sampled_columns <- sample(
      seq_len(ncol(expression_matrix)),
      size = ncol(expression_matrix),
      replace = TRUE
    )

    bootstrap_matrix <- expression_matrix[
      ,
      sampled_columns,
      drop = FALSE
    ]

    current <- conditioning_metrics_once(
      bootstrap_matrix,
      genes
    )

    current[
      ,
      bootstrap := bootstrap_index
    ]

    boot_rows[[bootstrap_index]] <- current
  }

  boot <- rbindlist(
    boot_rows,
    fill = TRUE
  )

  alpha <- (1 - confidence) / 2

  interval <- function(x) {
    finite <- x[is.finite(x)]

    if (length(finite) < MINIMUM_FINITE_BOOTSTRAPS) {
      return(
        c(
          low = NA_real_,
          high = NA_real_
        )
      )
    }

    stats::quantile(
      finite,
      probs = c(
        alpha,
        1 - alpha
      ),
      na.rm = TRUE,
      names = FALSE
    )
  }

  vif_mean_ci <- interval(boot$mean_vif)
  vif_max_ci <- interval(boot$max_vif)
  erank_ci <- interval(boot$effective_rank)
  normalised_erank_ci <- interval(boot$normalised_effective_rank)
  corr_ci <- interval(boot$mean_abs_correlation)

  summary <- copy(point)

  summary[
    ,
    `:=`(
      mean_vif_low = vif_mean_ci[[1]],
      mean_vif_high = vif_mean_ci[[2]],
      max_vif_low = vif_max_ci[[1]],
      max_vif_high = vif_max_ci[[2]],
      effective_rank_low = erank_ci[[1]],
      effective_rank_high = erank_ci[[2]],
      normalised_effective_rank_low = normalised_erank_ci[[1]],
      normalised_effective_rank_high = normalised_erank_ci[[2]],
      mean_abs_correlation_low = corr_ci[[1]],
      mean_abs_correlation_high = corr_ci[[2]],
      bootstraps_requested = n_boot,
      bootstraps_finite_erank = sum(is.finite(boot$effective_rank))
    )
  ]

  list(
    point = point,
    bootstrap = boot,
    summary = summary
  )
}


# ==============================================================================
# 10. STAGE 01 INPUT LOADERS
# ==============================================================================

load_full_development_bundle <- function(disease_id) {
  directory <- file.path(
    STAGE1_ROOT,
    "01_full_development",
    disease_id
  )

  matrix_path <- file.path(
    directory,
    paste0(
      disease_id,
      "_development_representation_matrices.rds"
    )
  )

  check_file_exists(
    matrix_path,
    paste0(disease_id, " Stage 01 full-development representations")
  )

  representations <- safe_read_rds(
    matrix_path,
    paste0(disease_id, " full-development representations")
  )

  metadata_candidates <- c(
    file.path(
      dirname(dirname(STAGE1_ROOT)),
      "data",
      "processed",
      disease_id
    ),
    file.path(
      "data",
      "processed",
      disease_id
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
    manifest_path,
    paste0(disease_id, " Stage 01 frozen manifest")
  )

  manifest <- as.data.frame(
    safe_read_csv(manifest_path),
    check.names = FALSE
  )

  development_manifest <- manifest[
    manifest$role == "development",
    ,
    drop = FALSE
  ]

  metadata_path <- file.path(
    "data",
    "processed",
    disease_id,
    "representation_specification",
    paste0(
      disease_id,
      "_canonical_metadata.rds"
    )
  )

  check_file_exists(
    metadata_path,
    paste0(disease_id, " canonical metadata")
  )

  metadata <- as.data.frame(
    safe_read_rds(
      metadata_path,
      paste0(disease_id, " canonical metadata")
    ),
    check.names = FALSE
  )

  metadata <- metadata[
    match(
      development_manifest$sample_id,
      metadata$sample_id
    ),
    ,
    drop = FALSE
  ]

  metadata$outcome <- factor(
    as.character(metadata$outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  list(
    directory = directory,
    representations = representations,
    metadata = metadata
  )
}

load_stage1_representation_pair <- function(directory,
                                            prefix_pattern = NULL) {
  development_path <- first_matching_file(
    directory,
    pattern = "_development_matrix\\.rds$",
    required = TRUE,
    label = "Stage 01 development representation matrix"
  )

  evaluation_path <- first_matching_file(
    directory,
    pattern = "_evaluation_matrix\\.rds$",
    required = FALSE,
    label = "Stage 01 evaluation representation matrix"
  )

  development_metadata_path <- first_matching_file(
    directory,
    pattern = "_development_metadata\\.csv$",
    required = TRUE,
    label = "Stage 01 development metadata"
  )

  evaluation_metadata_path <- first_matching_file(
    directory,
    pattern = "_evaluation_metadata\\.csv$",
    required = FALSE,
    label = "Stage 01 evaluation metadata"
  )

  development_matrix <- read_matrix_rds(
    development_path,
    "Stage 01 development representation matrix"
  )

  development_metadata <- read_metadata_csv(
    development_metadata_path,
    "Stage 01 development metadata"
  )

  aligned_development <- align_matrix_metadata(
    development_matrix,
    development_metadata,
    "Stage 01 development representation"
  )

  evaluation_matrix <- NULL
  evaluation_metadata <- NULL

  if (
    !is.na(evaluation_path) &&
      file.exists(evaluation_path)
  ) {
    evaluation_matrix <- read_matrix_rds(
      evaluation_path,
      "Stage 01 evaluation representation matrix"
    )

    if (
      is.na(evaluation_metadata_path) ||
        !file.exists(evaluation_metadata_path)
    ) {
      stop(
        "Evaluation matrix was found but its Stage 01 metadata file was absent in ",
        directory,
        call. = FALSE
      )
    }

    evaluation_metadata <- read_metadata_csv(
      evaluation_metadata_path,
      "Stage 01 evaluation metadata"
    )

    aligned_evaluation <- align_matrix_metadata(
      evaluation_matrix,
      evaluation_metadata,
      "Stage 01 evaluation representation"
    )

    evaluation_matrix <- aligned_evaluation$matrix
    evaluation_metadata <- aligned_evaluation$metadata

    assert_feature_alignment(
      aligned_development$matrix,
      evaluation_matrix,
      reference_label = "Stage 01 development representation",
      target_label = "Stage 01 evaluation representation"
    )

    assert_disjoint_samples(
      aligned_development$metadata$sample_id,
      evaluation_metadata$sample_id,
      development_label = "Stage 01 development",
      evaluation_label = "Stage 01 evaluation"
    )
  }

  list(
    development_matrix = aligned_development$matrix,
    development_metadata = aligned_development$metadata,
    evaluation_matrix = evaluation_matrix,
    evaluation_metadata = evaluation_metadata
  )
}


# ==============================================================================
# 10B. TRAINING-ONLY LOG2(CPM + 1) FOR ARCHITECTURE OUTER FOLDS
# ==============================================================================

fit_apply_log2_cpm_feature_stage <- function(train_counts,
                                             validation_counts = NULL,
                                             min_total_count = 10L,
                                             min_features = 50L) {
  training_total_counts <- rowSums(
    train_counts
  )

  keep_gene <- (
    is.finite(training_total_counts) &
      training_total_counts >= min_total_count
  )

  if (sum(keep_gene) < min_features) {
    stop(
      "Fewer than ",
      min_features,
      " genes survived the outer-training-only count filter.",
      call. = FALSE
    )
  }

  genes <- rownames(train_counts)[keep_gene]

  transform_one <- function(count_matrix) {
    count_matrix <- count_matrix[
      genes,
      ,
      drop = FALSE
    ]

    library_size <- colSums(
      count_matrix
    )

    if (
      any(!is.finite(library_size)) ||
        any(library_size <= 0)
    ) {
      stop(
        "Invalid library size during outer-fold log2(CPM + 1) transformation.",
        call. = FALSE
      )
    }

    cpm <- sweep(
      count_matrix,
      MARGIN = 2L,
      STATS = library_size,
      FUN = "/"
    ) * 1e6

    log2(
      cpm + 1
    )
  }

  list(
    train = transform_one(
      train_counts
    ),
    validation = if (!is.null(validation_counts)) {
      transform_one(
        validation_counts
      )
    } else {
      NULL
    },
    genes = genes
  )
}


# ==============================================================================
# 11. FULL-DEVELOPMENT FEATURE-SET CHARACTERISATION
# ==============================================================================

run_full_development_stage <- function(disease_id) {
  message_section(
    paste0(
      "Full-development feature-set characterisation | ",
      disease_id
    )
  )

  stage1 <- load_full_development_bundle(
    disease_id
  )

  output_dir <- file.path(
    STAGE2_ROOT,
    "01_full_development",
    disease_id
  )

  make_dir(output_dir)

  selected_set_rows <- list()
  jaccard_rows <- list()
  conditioning_rows <- list()
  cross_representation_rows <- list()

  selector_sets_by_representation <- list()

  representations <- stage1$representations
  metadata <- stage1$metadata

  representation_ids <- intersect(
    c(
      "RAW_counts",
      "VST_only",
      "VST_ComBat",
      "VST_ComBatRef",
      "VST_fSVA",
      "ComBatSeq_log1p"
    ),
    names(representations)
  )

  for (representation_index in seq_along(representation_ids)) {
    representation_id <- representation_ids[[representation_index]]

    message(
      "\n[representation] ",
      normalise_representation_label(representation_id)
    )

    matrix_current <- representations[[representation_id]]

    aligned <- align_matrix_metadata(
      matrix_current,
      metadata,
      paste0(
        disease_id,
        " ",
        representation_id
      )
    )

    matrix_current <- aligned$matrix
    metadata_current <- aligned$metadata

    representation_dir <- file.path(
      output_dir,
      safe_file_name(representation_id)
    )

    selector_outputs <- run_selector_collection(
      train_matrix = matrix_current,
      train_metadata = metadata_current,
      selector_ids = FULL_SELECTOR_IDS,
      top_n = DEFAULT_TOP_N,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "full_feature_characterisation",
        representation_id
      ),
      rf_trees = RF_IMPORTANCE_TREES_FULL,
      output_dir = file.path(
        representation_dir,
        "rankings"
      ),
      prefix = paste0(
        disease_id,
        "_",
        representation_id
      )
    )

    feature_sets <- lapply(
      selector_outputs,
      function(x) x$genes
    )

    selector_sets_by_representation[[representation_id]] <- feature_sets

    for (selector_id in names(feature_sets)) {
      genes <- feature_sets[[selector_id]]

      selected_set_rows[[length(selected_set_rows) + 1L]] <- data.table(
        disease = disease_id,
        representation = representation_id,
        representation_label = normalise_representation_label(representation_id),
        selector = selector_id,
        selector_label = normalise_selector_label(selector_id),
        top_n = DEFAULT_TOP_N,
        rank = seq_along(genes),
        gene = genes
      )

      conditioning <- bootstrap_conditioning_metrics(
        expression_matrix = matrix_current,
        genes = genes,
        n_boot = CONDITIONING_BOOTSTRAPS,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          representation_id,
          selector_id,
          "conditioning"
        )
      )

      conditioning_summary <- copy(
        conditioning$summary
      )

      conditioning_summary[
        ,
        `:=`(
          disease = disease_id,
          representation = representation_id,
          representation_label = normalise_representation_label(representation_id),
          selector = selector_id,
          selector_label = normalise_selector_label(selector_id),
          top_n = DEFAULT_TOP_N
        )
      ]

      conditioning_rows[[length(conditioning_rows) + 1L]] <- conditioning_summary

      bootstrap_path <- file.path(
        representation_dir,
        "conditioning_bootstraps",
        paste0(
          selector_id,
          "_conditioning_bootstrap.csv"
        )
      )

      safe_write_csv(
        conditioning$bootstrap,
        bootstrap_path
      )
    }

    matrix_ji <- jaccard_matrix(
      feature_sets
    )

    tidy_ji <- tidy_jaccard_matrix(
      matrix_object = matrix_ji,
      disease = disease_id,
      representation = representation_id,
      top_n = DEFAULT_TOP_N,
      context = "within_representation_across_selectors"
    )

    jaccard_rows[[representation_id]] <- tidy_ji

    safe_write_csv(
      as.data.table(
        matrix_ji,
        keep.rownames = "selector"
      ),
      file.path(
        representation_dir,
        paste0(
          disease_id,
          "_",
          representation_id,
          "_jaccard_matrix.csv"
        )
      )
    )

    plot_jaccard_heatmap(
      jaccard_dt = tidy_ji,
      output_path = file.path(
        representation_dir,
        paste0(
          disease_id,
          "_",
          representation_id,
          "_jaccard_heatmap.png"
        )
      ),
      title = paste0(
        disease_id,
        ": ",
        normalise_representation_label(representation_id),
        " feature-set concordance"
      )
    )
  }

  for (selector_id in FULL_SELECTOR_IDS) {
    available_representations <- names(
      selector_sets_by_representation
    )

    available_representations <- available_representations[
      vapply(
        selector_sets_by_representation[available_representations],
        function(x) selector_id %in% names(x),
        logical(1)
      )
    ]

    if (length(available_representations) < 2L) {
      next
    }

    for (i in seq_len(length(available_representations) - 1L)) {
      for (j in (i + 1L):length(available_representations)) {
        rep_1 <- available_representations[[i]]
        rep_2 <- available_representations[[j]]

        genes_1 <- selector_sets_by_representation[[rep_1]][[selector_id]]
        genes_2 <- selector_sets_by_representation[[rep_2]][[selector_id]]

        cross_representation_rows[[length(cross_representation_rows) + 1L]] <- data.table(
          disease = disease_id,
          selector = selector_id,
          selector_label = normalise_selector_label(selector_id),
          representation_1 = rep_1,
          representation_2 = rep_2,
          top_n = DEFAULT_TOP_N,
          shared_genes = length(intersect(genes_1, genes_2)),
          jaccard = jaccard_index(genes_1, genes_2)
        )
      }
    }
  }

  selected_set_table <- rbindlist(
    selected_set_rows,
    fill = TRUE
  )

  all_jaccard <- rbindlist(
    jaccard_rows,
    fill = TRUE
  )

  conditioning_table <- rbindlist(
    conditioning_rows,
    fill = TRUE
  )

  cross_representation_table <- rbindlist(
    cross_representation_rows,
    fill = TRUE
  )

  safe_write_csv(
    selected_set_table,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_full_development_feature_sets.csv"
      )
    )
  )

  safe_write_csv(
    all_jaccard,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_within_representation_jaccard.csv"
      )
    )
  )

  jaccard_summary <- summarise_jaccard(
    all_jaccard
  )

  if (is.list(jaccard_summary)) {
    safe_write_csv(
      jaccard_summary$representation,
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_jaccard_representation_summary.csv"
        )
      )
    )

    safe_write_csv(
      jaccard_summary$method,
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_jaccard_method_summary.csv"
        )
      )
    )
  }

  safe_write_csv(
    conditioning_table,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_feature_space_conditioning.csv"
      )
    )
  )

  safe_write_csv(
    cross_representation_table,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_same_selector_cross_representation_jaccard.csv"
      )
    )
  )

  p_conditioning <- ggplot(
    conditioning_table,
    aes(
      x = selector_label,
      y = normalised_effective_rank
    )
  ) +
    geom_point(size = 1.8) +
    geom_errorbar(
      aes(
        ymin = normalised_effective_rank_low,
        ymax = normalised_effective_rank_high
      ),
      width = 0.15
    ) +
    facet_wrap(
      ~ representation_label,
      scales = "free_x"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste0(
        disease_id,
        ": feature-space effective rank"
      ),
      x = "Feature-selection method",
      y = "Normalised effective rank"
    )

  ggsave(
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_normalised_effective_rank.png"
      )
    ),
    p_conditioning,
    width = 190,
    height = 125,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  invisible(
    list(
      feature_sets = selected_set_table,
      jaccard = all_jaccard,
      conditioning = conditioning_table,
      cross_representation = cross_representation_table
    )
  )
}


# ==============================================================================
# 12. EXPERIMENT 1: SAMPLE-SUPPORT FEATURE CHARACTERISATION
# ==============================================================================

run_experiment_1 <- function() {
  stage1_root <- file.path(
    STAGE1_ROOT,
    "02_experiment_1_sample_support"
  )

  stage2_root <- file.path(
    STAGE2_ROOT,
    "02_experiment_1_sample_support"
  )

  make_dir(stage2_root)

  experiment_diseases <- intersect(
    c("TB", "RA", "UC", "CD"),
    REQUESTED_DISEASES
  )

  for (disease_id in experiment_diseases) {
    message_section(
      paste0(
        "Experiment 1 feature characterisation | ",
        disease_id
      )
    )

    disease_stage1 <- file.path(
      stage1_root,
      disease_id,
      "01_representations"
    )

    disease_stage2 <- file.path(
      stage2_root,
      disease_id
    )

    make_dir(disease_stage2)

    representation_ids <- c(
      "VST_fSVA",
      "VST_ComBatRef",
      "ComBatSeq_log1p"
    )

    equal_size_sets <- list()
    all_conditioning <- list()

    for (representation_id in representation_ids) {
      pair_dir <- file.path(
        disease_stage1,
        representation_id
      )

      pair <- load_stage1_representation_pair(
        pair_dir
      )

      representation_output <- file.path(
        disease_stage2,
        representation_id
      )

      selector_outputs <- run_selector_collection(
        train_matrix = pair$development_matrix,
        train_metadata = pair$development_metadata,
        selector_ids = EXPERIMENT1_SELECTOR_IDS,
        top_n = 50L,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment1",
          representation_id
        ),
        rf_trees = RF_IMPORTANCE_TREES_FULL,
        output_dir = file.path(
          representation_output,
          "rankings"
        ),
        prefix = paste0(
          disease_id,
          "_equal_size_",
          representation_id
        )
      )

      feature_sets <- lapply(
        selector_outputs,
        function(x) x$genes
      )

      equal_size_sets[[representation_id]] <- feature_sets

      for (selector_id in names(selector_outputs)) {
        result <- selector_outputs[[selector_id]]

        selector_dir <- file.path(
          representation_output,
          "selected_matrices",
          selector_id
        )

        save_selected_matrix_pair(
          train_matrix = pair$development_matrix,
          evaluation_matrix = pair$evaluation_matrix,
          genes = result$genes,
          output_dir = selector_dir,
          prefix = paste0(
            disease_id,
            "_equal_size_",
            representation_id,
            "_",
            selector_id,
            "_top50"
          )
        )

        conditioning <- bootstrap_conditioning_metrics(
          expression_matrix = pair$development_matrix,
          genes = result$genes,
          n_boot = CONDITIONING_BOOTSTRAPS,
          seed = make_trace_subseed(
            GLOBAL_SEED,
            disease_id,
            "experiment1",
            representation_id,
            selector_id,
            "conditioning"
          )
        )

        current_conditioning <- copy(
          conditioning$summary
        )

        current_conditioning[
          ,
          `:=`(
            disease = disease_id,
            representation = representation_id,
            selector = selector_id,
            top_n = 50L,
            context = "equal_size_development"
          )
        ]

        all_conditioning[[length(all_conditioning) + 1L]] <- current_conditioning

        safe_write_csv(
          conditioning$bootstrap,
          file.path(
            representation_output,
            "conditioning",
            paste0(
              selector_id,
              "_bootstrap.csv"
            )
          )
        )
      }

      ji <- jaccard_matrix(
        feature_sets
      )

      ji_tidy <- tidy_jaccard_matrix(
        matrix_object = ji,
        disease = disease_id,
        representation = representation_id,
        top_n = 50L,
        context = "experiment1_equal_size"
      )

      safe_write_csv(
        ji_tidy,
        file.path(
          representation_output,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_equal_size_jaccard.csv"
          )
        )
      )

      plot_jaccard_heatmap(
        ji_tidy,
        file.path(
          representation_output,
          paste0(
            disease_id,
            "_",
            representation_id,
            "_equal_size_jaccard.png"
          )
        ),
        paste0(
          disease_id,
          ": equal-size feature-set concordance - ",
          normalise_representation_label(representation_id)
        )
      )
    }

    conditioning_dt <- rbindlist(
      all_conditioning,
      fill = TRUE
    )

    safe_write_csv(
      conditioning_dt,
      file.path(
        disease_stage2,
        paste0(
          disease_id,
          "_equal_size_conditioning.csv"
        )
      )
    )

    full_stage2_dir <- file.path(
      STAGE2_ROOT,
      "01_full_development",
      disease_id
    )

    full_feature_path <- file.path(
      full_stage2_dir,
      paste0(
        disease_id,
        "_full_development_feature_sets.csv"
      )
    )

    if (file.exists(full_feature_path)) {
      full_features <- as.data.table(
        safe_read_csv(
          full_feature_path
        )
      )

      retention_rows <- list()

      for (representation_id in names(equal_size_sets)) {
        for (selector_id in names(equal_size_sets[[representation_id]])) {
          full_genes <- full_features[
            representation == representation_id &
              selector == selector_id,
            gene
          ]

          equal_genes <- equal_size_sets[[representation_id]][[selector_id]]

          if (length(full_genes) == 0L) {
            next
          }

          retention_rows[[length(retention_rows) + 1L]] <- data.table(
            disease = disease_id,
            representation = representation_id,
            selector = selector_id,
            full_top_n = length(full_genes),
            equal_size_top_n = length(equal_genes),
            shared_genes = length(
              intersect(
                full_genes,
                equal_genes
              )
            ),
            jaccard = jaccard_index(
              full_genes,
              equal_genes
            )
          )
        }
      }

      retention_dt <- rbindlist(
        retention_rows,
        fill = TRUE
      )

      safe_write_csv(
        retention_dt,
        file.path(
          disease_stage2,
          paste0(
            disease_id,
            "_full_vs_equal_size_feature_retention.csv"
          )
        )
      )
    }
  }

  invisible(TRUE)
}


# ==============================================================================
# 13. EXPERIMENT 2: DISEASE-COMPOSITION FEATURE CHARACTERISATION
# ==============================================================================

run_experiment_2 <- function() {
  if (!all(c("TB", "RA", "UC", "CD", "T1D") %in% REQUESTED_DISEASES)) {
    warning(
      "Experiment 2 requires --diseases=TB,RA,UC,CD,T1D. Skipping.",
      call. = FALSE
    )

    return(invisible(NULL))
  }

  stage1_root <- file.path(
    STAGE1_ROOT,
    "03_experiment_2_disease_composition"
  )

  stage2_root <- file.path(
    STAGE2_ROOT,
    "03_experiment_2_disease_composition"
  )

  make_dir(stage2_root)

  message_section(
    "Experiment 2 feature characterisation | disease composition"
  )

  model_directories <- list(
    Model_A = file.path(
      stage1_root,
      "Model_A_TB_naive"
    ),
    Model_B = file.path(
      stage1_root,
      "Model_B_ALL5"
    )
  )

  model_feature_sets <- list()

  for (model_name in names(model_directories)) {
    pair <- load_stage1_representation_pair(
      model_directories[[model_name]]
    )

    model_output <- file.path(
      stage2_root,
      model_name
    )

    selector_outputs <- run_selector_collection(
      train_matrix = pair$development_matrix,
      train_metadata = pair$development_metadata,
      selector_ids = EXPERIMENT2_MODEL_SELECTOR_IDS,
      top_n = 50L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "experiment2",
        model_name
      ),
      rf_trees = RF_IMPORTANCE_TREES_FULL,
      output_dir = file.path(
        model_output,
        "rankings"
      ),
      prefix = paste0(
        model_name,
        "_VST_fSVA"
      )
    )

    model_feature_sets[[model_name]] <- lapply(
      selector_outputs,
      function(x) x$genes
    )

    for (selector_id in names(selector_outputs)) {
      save_selected_matrix_pair(
        train_matrix = pair$development_matrix,
        evaluation_matrix = pair$evaluation_matrix,
        genes = selector_outputs[[selector_id]]$genes,
        output_dir = file.path(
          model_output,
          "selected_matrices",
          selector_id
        ),
        prefix = paste0(
          model_name,
          "_",
          selector_id,
          "_top50"
        )
      )
    }
  }

  comparison_rows <- list()

  for (selector_id in EXPERIMENT2_MODEL_SELECTOR_IDS) {
    genes_a <- model_feature_sets[["Model_A"]][[selector_id]]
    genes_b <- model_feature_sets[["Model_B"]][[selector_id]]

    comparison_rows[[selector_id]] <- data.table(
      selector = selector_id,
      selector_label = normalise_selector_label(selector_id),
      model_a_top_n = length(genes_a),
      model_b_top_n = length(genes_b),
      shared_genes = length(intersect(genes_a, genes_b)),
      jaccard = jaccard_index(genes_a, genes_b)
    )
  }

  safe_write_csv(
    rbindlist(
      comparison_rows,
      fill = TRUE
    ),
    file.path(
      stage2_root,
      "Model_A_vs_Model_B_feature_concordance.csv"
    )
  )

  lodo_root <- file.path(
    stage1_root,
    "LODO"
  )

  target_directories <- sort(
    list.dirs(
      lodo_root,
      full.names = TRUE,
      recursive = FALSE
    )
  )

  lodo_rows <- list()

  for (target_directory in target_directories) {
    target_name <- sub(
      "^target_",
      "",
      basename(target_directory)
    )

    pair <- load_stage1_representation_pair(
      target_directory
    )

    result <- select_ranked_features(
      train_matrix = pair$development_matrix,
      train_metadata = pair$development_metadata,
      filter_id = EXPERIMENT2_LOCKED_SELECTOR_ID,
      top_n = 50L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "experiment2",
        "LODO",
        target_name,
        EXPERIMENT2_LOCKED_SELECTOR_ID
      )
    )

    target_output <- file.path(
      stage2_root,
      "LODO",
      paste0(
        "target_",
        target_name
      )
    )

    make_dir(target_output)

    safe_write_csv(
      result$scores,
      file.path(
        target_output,
        paste0(
          "LODO_",
          target_name,
          "_mRMR_full_ranking.csv"
        )
      )
    )

    save_gene_vector(
      result$genes,
      file.path(
        target_output,
        paste0(
          "LODO_",
          target_name,
          "_mRMR_top50_genes.rds"
        )
      )
    )

    save_selected_matrix_pair(
      train_matrix = pair$development_matrix,
      evaluation_matrix = pair$evaluation_matrix,
      genes = result$genes,
      output_dir = target_output,
      prefix = paste0(
        "LODO_",
        target_name,
        "_mRMR_top50"
      )
    )

    lodo_rows[[target_name]] <- data.table(
      target_disease = target_name,
      selector = "mrmr",
      top_n = 50L,
      source_n = ncol(pair$development_matrix),
      target_n = ncol(pair$evaluation_matrix),
      shared_feature_space_before_selection = nrow(pair$development_matrix)
    )
  }

  safe_write_csv(
    rbindlist(
      lodo_rows,
      fill = TRUE
    ),
    file.path(
      stage2_root,
      "LODO_feature_selection_summary.csv"
    )
  )

  invisible(TRUE)
}


# ==============================================================================
# 14. EXPERIMENT 3: COHORT-COMPOSITION / LOCO FEATURE CHARACTERISATION
# ==============================================================================

run_rf_candidate_feature_counts <- function(pair,
                                            output_dir,
                                            prefix,
                                            seed,
                                            rf_trees) {
  max_top_n <- max(
    EXPERIMENT3_FEATURE_COUNTS
  )

  ranking <- select_ranked_features(
    train_matrix = pair$development_matrix,
    train_metadata = pair$development_metadata,
    filter_id = "importance",
    top_n = max_top_n,
    seed = seed,
    rf_trees = rf_trees
  )

  make_dir(output_dir)

  safe_write_csv(
    ranking$scores,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_RF_importance_full_ranking.csv"
      )
    )
  )

  candidate_rows <- list()

  for (feature_count in EXPERIMENT3_FEATURE_COUNTS) {
    genes <- ranking$scores[
      seq_len(feature_count),
      gene
    ]

    save_gene_vector(
      genes,
      file.path(
        output_dir,
        paste0(
          prefix,
          "_RF_importance_top",
          feature_count,
          "_genes.rds"
        )
      )
    )

    save_selected_matrix_pair(
      train_matrix = pair$development_matrix,
      evaluation_matrix = pair$evaluation_matrix,
      genes = genes,
      output_dir = file.path(
        output_dir,
        paste0(
          "top",
          feature_count
        )
      ),
      prefix = paste0(
        prefix,
        "_RF_importance_top",
        feature_count
      )
    )

    candidate_rows[[as.character(feature_count)]] <- data.table(
      top_n = feature_count,
      genes_available = nrow(pair$development_matrix),
      development_n = ncol(pair$development_matrix),
      evaluation_n = if (!is.null(pair$evaluation_matrix)) ncol(pair$evaluation_matrix) else 0L
    )
  }

  list(
    ranking = ranking,
    summary = rbindlist(
      candidate_rows,
      fill = TRUE
    )
  )
}

run_experiment_3 <- function() {
  stage1_root <- file.path(
    STAGE1_ROOT,
    "04_experiment_3_cohort_composition"
  )

  stage2_root <- file.path(
    STAGE2_ROOT,
    "04_experiment_3_cohort_composition"
  )

  make_dir(stage2_root)

  for (disease_id in REQUESTED_DISEASES) {
    disease_stage1 <- file.path(
      stage1_root,
      disease_id
    )

    if (!dir.exists(disease_stage1)) {
      warning(
        "Experiment 3 Stage 01 directory absent for ",
        disease_id,
        ". Skipping.",
        call. = FALSE
      )

      next
    }

    disease_stage2 <- file.path(
      stage2_root,
      disease_id
    )

    make_dir(disease_stage2)

    message_section(
      paste0(
        "Experiment 3 feature characterisation | ",
        disease_id
      )
    )

    representation_ids <- c(
      "VST_only",
      "log2_CPM_plus_1"
    )

    pooled_rows <- list()

    for (representation_id in representation_ids) {
      pair_dir <- file.path(
        disease_stage1,
        "frozen_heldout",
        representation_id
      )

      pair <- load_stage1_representation_pair(
        pair_dir
      )

      result <- run_rf_candidate_feature_counts(
        pair = pair,
        output_dir = file.path(
          disease_stage2,
          "frozen_heldout",
          representation_id
        ),
        prefix = paste0(
          disease_id,
          "_",
          representation_id
        ),
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment3",
          "frozen_heldout",
          representation_id
        ),
        rf_trees = RF_IMPORTANCE_TREES_DISEASE_SPECIFIC
      )

      current_summary <- copy(
        result$summary
      )

      current_summary[
        ,
        `:=`(
          disease = disease_id,
          representation = representation_id,
          context = "frozen_heldout"
        )
      ]

      pooled_rows[[representation_id]] <- current_summary
    }

    safe_write_csv(
      rbindlist(
        pooled_rows,
        fill = TRUE
      ),
      file.path(
        disease_stage2,
        paste0(
          disease_id,
          "_frozen_heldout_feature_count_summary.csv"
        )
      )
    )

    loco_root <- file.path(
      disease_stage1,
      "LOCO"
    )

    target_cohort_directories <- sort(
      list.dirs(
        loco_root,
        full.names = TRUE,
        recursive = FALSE
      )
    )

    loco_rows <- list()

    for (target_directory in target_cohort_directories) {
      target_cohort <- basename(
        target_directory
      )

      for (representation_id in representation_ids) {
        pair_dir <- file.path(
          target_directory,
          representation_id
        )

        if (!dir.exists(pair_dir)) {
          next
        }

        pair <- load_stage1_representation_pair(
          pair_dir
        )

        result <- run_rf_candidate_feature_counts(
          pair = pair,
          output_dir = file.path(
            disease_stage2,
            "LOCO",
            target_cohort,
            representation_id
          ),
          prefix = paste0(
            disease_id,
            "_LOCO_",
            target_cohort,
            "_",
            representation_id
          ),
          seed = make_trace_subseed(
            GLOBAL_SEED,
            disease_id,
            "experiment3",
            "LOCO",
            target_cohort,
            representation_id
          ),
          rf_trees = RF_IMPORTANCE_TREES_DISEASE_SPECIFIC
        )

        current_summary <- copy(
          result$summary
        )

        current_summary[
          ,
          `:=`(
            disease = disease_id,
            target_cohort = target_cohort,
            representation = representation_id,
            context = "LOCO"
          )
        ]

        loco_rows[[paste(target_cohort, representation_id, sep = "__")]] <- current_summary
      }
    }

    safe_write_csv(
      rbindlist(
        loco_rows,
        fill = TRUE
      ),
      file.path(
        disease_stage2,
        paste0(
          disease_id,
          "_LOCO_feature_count_summary.csv"
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 15. EXPERIMENT 4: MATCHED ARCHITECTURE FEATURE CHARACTERISATION
# ==============================================================================

run_experiment_4 <- function() {
  stage1_root <- file.path(
    STAGE1_ROOT,
    "05_experiment_4_model_architecture"
  )

  stage2_root <- file.path(
    STAGE2_ROOT,
    "05_experiment_4_model_architecture"
  )

  make_dir(stage2_root)

  for (disease_id in REQUESTED_DISEASES) {
    disease_stage1 <- file.path(
      stage1_root,
      disease_id
    )

    if (!dir.exists(disease_stage1)) {
      warning(
        "Experiment 4 Stage 01 directory absent for ",
        disease_id,
        ". Skipping.",
        call. = FALSE
      )

      next
    }

    disease_stage2 <- file.path(
      stage2_root,
      disease_id
    )

    make_dir(disease_stage2)

    message_section(
      paste0(
        "Experiment 4 matched architecture feature selection | ",
        disease_id
      )
    )

    frozen_pair <- load_stage1_representation_pair(
      file.path(
        disease_stage1,
        "frozen_diagnostic_test"
      )
    )

    final_ranking <- select_ranked_features(
      train_matrix = frozen_pair$development_matrix,
      train_metadata = frozen_pair$development_metadata,
      filter_id = EXPERIMENT4_SELECTOR_ID,
      top_n = EXPERIMENT4_TOP_N,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "experiment4",
        "frozen_diagnostic_test"
      ),
      rf_trees = RF_IMPORTANCE_TREES_ARCHITECTURE
    )

    final_dir <- file.path(
      disease_stage2,
      "frozen_diagnostic_test"
    )

    make_dir(final_dir)

    safe_write_csv(
      final_ranking$scores,
      file.path(
        final_dir,
        paste0(
          disease_id,
          "_architecture_RF_importance_full_ranking.csv"
        )
      )
    )

    save_gene_vector(
      final_ranking$genes,
      file.path(
        final_dir,
        paste0(
          disease_id,
          "_architecture_top",
          EXPERIMENT4_TOP_N,
          "_genes.rds"
        )
      )
    )

    save_selected_matrix_pair(
      train_matrix = frozen_pair$development_matrix,
      evaluation_matrix = frozen_pair$evaluation_matrix,
      genes = final_ranking$genes,
      output_dir = final_dir,
      prefix = paste0(
        disease_id,
        "_architecture_RF_importance_top",
        EXPERIMENT4_TOP_N
      )
    )

    matched_architecture_audit <- data.table(
      disease = disease_id,
      selector = "RF permutation importance",
      top_n = EXPERIMENT4_TOP_N,
      rf_importance_trees = RF_IMPORTANCE_TREES_ARCHITECTURE,
      development_n = ncol(frozen_pair$development_matrix),
      frozen_test_n = ncol(frozen_pair$evaluation_matrix),
      LightGBM_gene_vector = basename(
        file.path(
          final_dir,
          paste0(
            disease_id,
            "_architecture_top",
            EXPERIMENT4_TOP_N,
            "_genes.rds"
          )
        )
      ),
      TabPFN_gene_vector = basename(
        file.path(
          final_dir,
          paste0(
            disease_id,
            "_architecture_top",
            EXPERIMENT4_TOP_N,
            "_genes.rds"
          )
        )
      ),
      identical_gene_vector = TRUE
    )

    safe_write_csv(
      matched_architecture_audit,
      file.path(
        final_dir,
        paste0(
          disease_id,
          "_matched_architecture_feature_audit.csv"
        )
      )
    )

    fold_map_path <- file.path(
      disease_stage1,
      paste0(
        disease_id,
        "_architecture_outer_fold_map.csv"
      )
    )

    if (file.exists(fold_map_path)) {
      fold_map <- as.data.frame(
        safe_read_csv(
          fold_map_path
        ),
        check.names = FALSE
      )

      check_required_columns(
        fold_map,
        c(
          "sample_id",
          "fold"
        ),
        paste0(
          disease_id,
          " architecture outer fold map"
        )
      )

      outer_fold_rows <- list()
      fold_feature_sets <- list()

      for (fold_id in sort(unique(fold_map$fold))) {
        validation_ids <- fold_map$sample_id[
          fold_map$fold == fold_id
        ]

        training_ids <- setdiff(
          fold_map$sample_id,
          validation_ids
        )

        assert_disjoint_samples(
          training_ids,
          validation_ids,
          development_label = paste0(
            "architecture outer-fold ",
            fold_id,
            " training"
          ),
          evaluation_label = paste0(
            "architecture outer-fold ",
            fold_id,
            " validation"
          )
        )

        canonical_counts_path <- file.path(
          "data",
          "processed",
          disease_id,
          "representation_specification",
          paste0(
            disease_id,
            "_canonical_raw_counts.rds"
          )
        )

        check_file_exists(
          canonical_counts_path,
          paste0(
            disease_id,
            " canonical raw counts for architecture outer folds"
          )
        )

        canonical_counts <- read_matrix_rds(
          canonical_counts_path,
          paste0(
            disease_id,
            " canonical raw counts"
          )
        )

        fold_train_counts <- canonical_counts[
          ,
          training_ids,
          drop = FALSE
        ]

        fold_validation_counts <- canonical_counts[
          ,
          validation_ids,
          drop = FALSE
        ]

        fold_representation <- fit_apply_log2_cpm_feature_stage(
          train_counts = fold_train_counts,
          validation_counts = fold_validation_counts,
          min_total_count = as.integer(
            get_config_value(
              config,
              "analysis.min_gene_count",
              10
            )
          ),
          min_features = EXPERIMENT4_TOP_N
        )

        fold_train_matrix <- fold_representation$train
        fold_validation_matrix <- fold_representation$validation

        fold_train_metadata <- frozen_pair$development_metadata[
          match(
            training_ids,
            frozen_pair$development_metadata$sample_id
          ),
          ,
          drop = FALSE
        ]

        fold_ranking <- select_ranked_features(
          train_matrix = fold_train_matrix,
          train_metadata = fold_train_metadata,
          filter_id = EXPERIMENT4_SELECTOR_ID,
          top_n = EXPERIMENT4_TOP_N,
          seed = make_trace_subseed(
            GLOBAL_SEED,
            disease_id,
            "experiment4",
            "outer_fold",
            fold_id
          ),
          rf_trees = RF_IMPORTANCE_TREES_ARCHITECTURE
        )

        fold_dir <- file.path(
          disease_stage2,
          "outer_cv",
          paste0(
            "fold_",
            fold_id
          )
        )

        make_dir(fold_dir)

        safe_write_csv(
          fold_ranking$scores,
          file.path(
            fold_dir,
            paste0(
              disease_id,
              "_outer_fold_",
              fold_id,
              "_RF_importance_full_ranking.csv"
            )
          )
        )

        save_gene_vector(
          fold_ranking$genes,
          file.path(
            fold_dir,
            paste0(
              disease_id,
              "_outer_fold_",
              fold_id,
              "_top",
              EXPERIMENT4_TOP_N,
              "_genes.rds"
            )
          )
        )

        save_selected_matrix_pair(
          train_matrix = fold_train_matrix,
          evaluation_matrix = fold_validation_matrix,
          genes = fold_ranking$genes,
          output_dir = fold_dir,
          prefix = paste0(
            disease_id,
            "_outer_fold_",
            fold_id,
            "_RF_importance_top",
            EXPERIMENT4_TOP_N
          )
        )

        fold_feature_sets[[as.character(fold_id)]] <- fold_ranking$genes

        outer_fold_rows[[as.character(fold_id)]] <- data.table(
          disease = disease_id,
          fold = fold_id,
          train_n = length(training_ids),
          validation_n = length(validation_ids),
          top_n = EXPERIMENT4_TOP_N,
          selector = "importance",
          rf_importance_trees = RF_IMPORTANCE_TREES_ARCHITECTURE
        )
      }

      safe_write_csv(
        rbindlist(
          outer_fold_rows,
          fill = TRUE
        ),
        file.path(
          disease_stage2,
          paste0(
            disease_id,
            "_architecture_outer_fold_feature_summary.csv"
          )
        )
      )

      if (length(fold_feature_sets) >= 2L) {
        fold_ji <- jaccard_matrix(
          fold_feature_sets
        )

        fold_ji_dt <- tidy_jaccard_matrix(
          matrix_object = fold_ji,
          disease = disease_id,
          representation = "log2_CPM_plus_1",
          top_n = EXPERIMENT4_TOP_N,
          context = "architecture_outer_fold_feature_stability"
        )

        safe_write_csv(
          fold_ji_dt,
          file.path(
            disease_stage2,
            paste0(
              disease_id,
              "_architecture_outer_fold_jaccard.csv"
            )
          )
        )
      }
    }

    loco_root <- file.path(
      disease_stage1,
      "LOCO"
    )

    if (dir.exists(loco_root)) {
      target_cohort_directories <- sort(
        list.dirs(
          loco_root,
          full.names = TRUE,
          recursive = FALSE
        )
      )

      loco_feature_sets <- list()
      loco_rows <- list()

      for (target_directory in target_cohort_directories) {
        target_cohort <- basename(
          target_directory
        )

        pair <- load_stage1_representation_pair(
          target_directory
        )

        ranking <- select_ranked_features(
          train_matrix = pair$development_matrix,
          train_metadata = pair$development_metadata,
          filter_id = EXPERIMENT4_SELECTOR_ID,
          top_n = EXPERIMENT4_TOP_N,
          seed = make_trace_subseed(
            GLOBAL_SEED,
            disease_id,
            "experiment4",
            "LOCO",
            target_cohort
          ),
          rf_trees = RF_IMPORTANCE_TREES_ARCHITECTURE
        )

        target_output <- file.path(
          disease_stage2,
          "LOCO",
          target_cohort
        )

        make_dir(target_output)

        safe_write_csv(
          ranking$scores,
          file.path(
            target_output,
            paste0(
              disease_id,
              "_architecture_LOCO_",
              target_cohort,
              "_RF_importance_full_ranking.csv"
            )
          )
        )

        save_gene_vector(
          ranking$genes,
          file.path(
            target_output,
            paste0(
              disease_id,
              "_architecture_LOCO_",
              target_cohort,
              "_top",
              EXPERIMENT4_TOP_N,
              "_genes.rds"
            )
          )
        )

        save_selected_matrix_pair(
          train_matrix = pair$development_matrix,
          evaluation_matrix = pair$evaluation_matrix,
          genes = ranking$genes,
          output_dir = target_output,
          prefix = paste0(
            disease_id,
            "_architecture_LOCO_",
            target_cohort,
            "_RF_importance_top",
            EXPERIMENT4_TOP_N
          )
        )

        loco_feature_sets[[target_cohort]] <- ranking$genes

        loco_rows[[target_cohort]] <- data.table(
          disease = disease_id,
          target_cohort = target_cohort,
          source_n = ncol(pair$development_matrix),
          target_n = ncol(pair$evaluation_matrix),
          top_n = EXPERIMENT4_TOP_N,
          selector = "importance",
          rf_importance_trees = RF_IMPORTANCE_TREES_ARCHITECTURE
        )
      }

      safe_write_csv(
        rbindlist(
          loco_rows,
          fill = TRUE
        ),
        file.path(
          disease_stage2,
          paste0(
            disease_id,
            "_architecture_LOCO_feature_summary.csv"
          )
        )
      )

      if (length(loco_feature_sets) >= 2L) {
        loco_ji <- jaccard_matrix(
          loco_feature_sets
        )

        loco_ji_dt <- tidy_jaccard_matrix(
          matrix_object = loco_ji,
          disease = disease_id,
          representation = "log2_CPM_plus_1",
          top_n = EXPERIMENT4_TOP_N,
          context = "architecture_LOCO_feature_stability"
        )

        safe_write_csv(
          loco_ji_dt,
          file.path(
            disease_stage2,
            paste0(
              disease_id,
              "_architecture_LOCO_jaccard.csv"
            )
          )
        )
      }
    }
  }

  invisible(TRUE)
}


# ==============================================================================
# 16. CROSS-DISEASE / CROSS-CONTEXT AUDIT TABLES
# ==============================================================================

write_stage2_method_manifest <- function() {
  selector_manifest <- data.table(
    selector_id = FULL_SELECTOR_IDS,
    selector_label = unname(FULL_SELECTOR_LABELS[FULL_SELECTOR_IDS]),
    full_development = TRUE,
    experiment_1 = FULL_SELECTOR_IDS %in% EXPERIMENT1_SELECTOR_IDS,
    experiment_2_model_ab = FULL_SELECTOR_IDS %in% EXPERIMENT2_MODEL_SELECTOR_IDS,
    experiment_2_locked_lodo = FULL_SELECTOR_IDS == EXPERIMENT2_LOCKED_SELECTOR_ID,
    experiment_3 = FULL_SELECTOR_IDS == EXPERIMENT3_SELECTOR_ID,
    experiment_4 = FULL_SELECTOR_IDS == EXPERIMENT4_SELECTOR_ID
  )

  safe_write_csv(
    selector_manifest,
    file.path(
      AUDIT_ROOT,
      "TRACE_Omics_stage02_selector_manifest.csv"
    )
  )

  feature_count_manifest <- data.table(
    context = c(
      "full_development",
      "experiment_1_equal_size",
      "experiment_2_model_A_B",
      "experiment_2_LODO",
      "experiment_3_LOCO_candidates",
      "experiment_4_architecture"
    ),
    feature_counts = c(
      as.character(DEFAULT_TOP_N),
      "50",
      "50",
      "50",
      paste(EXPERIMENT3_FEATURE_COUNTS, collapse = "/"),
      as.character(EXPERIMENT4_TOP_N)
    )
  )

  safe_write_csv(
    feature_count_manifest,
    file.path(
      AUDIT_ROOT,
      "TRACE_Omics_stage02_feature_count_manifest.csv"
    )
  )
}

write_stage2_method_manifest()


# ==============================================================================
# 17. RUN REQUESTED WORKFLOWS
# ==============================================================================

if (run_requested("full")) {
  for (disease_id in REQUESTED_DISEASES) {
    run_full_development_stage(
      disease_id
    )
  }
}

if (run_requested("experiment1")) {
  run_experiment_1()
}

if (run_requested("experiment2")) {
  run_experiment_2()
}

if (run_requested("experiment3")) {
  run_experiment_3()
}

if (run_requested("experiment4")) {
  run_experiment_4()
}


# ==============================================================================
# 18. FINAL REPRODUCIBILITY RECORDS
# ==============================================================================

config_snapshot_path <- file.path(
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
  config_snapshot_path,
  overwrite = TRUE
)

write_session_info(
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage02_feature_set_characterisation_sessionInfo.txt"
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
  full_development_top_n = DEFAULT_TOP_N,
  conditioning_bootstraps = CONDITIONING_BOOTSTRAPS,
  rf_importance_trees_full = RF_IMPORTANCE_TREES_FULL,
  rf_importance_trees_disease_specific = RF_IMPORTANCE_TREES_DISEASE_SPECIFIC,
  rf_importance_trees_architecture = RF_IMPORTANCE_TREES_ARCHITECTURE,
  strict_audits = STRICT_AUDITS
)

safe_write_csv(
  run_manifest,
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage02_run_manifest.csv"
  )
)

if (CONDITIONING_BOOTSTRAPS != 200L) {
  warning(
    "The source Stage 02 workflow supplied for this rebuild used 200 conditioning ",
    "bootstrap iterations. The current configuration resolves to ",
    CONDITIONING_BOOTSTRAPS,
    ". Set feature_set_characterisation.conditioning.bootstrap_iterations: 200 ",
    "for exact reproduction of that analysis.",
    call. = FALSE
  )
}

cat(
  "\n============================================================\n",
  "TRACE-OMICS STAGE 02 COMPLETE\n",
  "Feature-set characterisation finished successfully.\n",
  "Output root: ",
  STAGE2_ROOT,
  "\n============================================================\n",
  sep = ""
)

# ==============================================================================
# End of scripts/02_feature_set_characterisation/02_feature_set_characterisation.R
# ==============================================================================
