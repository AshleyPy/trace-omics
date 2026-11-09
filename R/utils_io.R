# ==============================================================================
# TRACE-Omics
# Transcriptomic Robustness Assessment and Cross-cohort Evaluation for Omics
# Core utility functions for configuration, validation, task handling and audit
#
# File:
#   R/utils_io.R
#
# Purpose:
#   Reusable helpers shared across the four TRACE-Omics analytical stages:
#   representation specification, feature-set characterisation, model development
#   and optimisation, and generalisation and reliability assessment. The helpers
#   keep sample identity, cohort structure, feature order, configuration state and
#   prediction outputs explicit so that the controlled perturbation experiments can
#   be reproduced and audited consistently.
# ==============================================================================

# ==============================================================================
# 1. Package helpers
# ==============================================================================

require_packages <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    msg <- paste0(
      "The following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nInstall them before running the workflow."
    )

    if (isTRUE(bioc)) {
      msg <- paste0(
        msg,
        "\nFor Bioconductor packages, use:\n",
        "  BiocManager::install(c(",
        paste(sprintf('"%s"', missing), collapse = ", "),
        "))"
      )
    } else {
      msg <- paste0(
        msg,
        "\nFor CRAN packages, use:\n",
        "  install.packages(c(",
        paste(sprintf('"%s"', missing), collapse = ", "),
        "))"
      )
    }

    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}

quiet_library <- function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

# ==============================================================================
# 2. Path helpers
# ==============================================================================

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

check_file_exists <- function(path, label = "File") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }

  invisible(path)
}

check_dir_exists <- function(path, label = "Directory") {
  if (!dir.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }

  invisible(path)
}

safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

path_here <- function(...) {
  file.path(...)
}

write_session_info <- function(path) {
  make_dir(dirname(path))

  writeLines(
    capture.output(sessionInfo()),
    con = path
  )

  invisible(path)
}

safe_save_rds <- function(object, path, compress = "gzip") {
  make_dir(dirname(path))
  saveRDS(object, file = path, compress = compress)
  invisible(path)
}

safe_read_rds <- function(path, label = "RDS file") {
  check_file_exists(path, label)
  readRDS(path)
}

safe_write_csv <- function(x, path) {
  make_dir(dirname(path))

  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(as.data.frame(x), path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }

  invisible(path)
}

safe_read_csv <- function(path, ...) {
  check_file_exists(path, "CSV file")

  if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, ...)
  } else {
    utils::read.csv(path, check.names = FALSE, ...)
  }
}

# ==============================================================================
# 3. Configuration helpers
# ==============================================================================

read_trace_config <- function(path = "config/example_config.yml") {
  check_file_exists(path, "Configuration file")

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "Package 'yaml' is required to read configuration files.\n",
      "Install it with: install.packages(\"yaml\")",
      call. = FALSE
    )
  }

  yaml::read_yaml(path)
}

get_config_value <- function(config, path, default = NULL) {
  parts <- strsplit(path, "\\.")[[1]]

  value <- config

  for (part in parts) {
    if (is.null(value[[part]])) {
      return(default)
    }

    value <- value[[part]]
  }

  value
}

print_config_summary <- function(config) {
  disease_id <- get_config_value(config, "analysis.disease_id", "unknown")
  seed <- get_config_value(config, "project.seed", NA)
  positive_class <- get_config_value(config, "analysis.positive_class", "case")
  negative_class <- get_config_value(config, "analysis.negative_class", "control")

  message("TRACE-Omics configuration")
  message("  Disease ID:      ", disease_id)
  message("  Seed:            ", seed)
  message("  Positive class:  ", positive_class)
  message("  Negative class:  ", negative_class)

  invisible(config)
}

trace_stages <- function() {
  c(
    "Representation specification",
    "Feature-set characterisation",
    "Model development and optimisation",
    "Generalisation and reliability assessment"
  )
}

trace_experiments <- function() {
  c(
    "Sample support",
    "Disease composition",
    "Cohort composition",
    "Model architecture"
  )
}

print_trace_overview <- function() {
  message("TRACE-Omics analytical stages")
  for (i in seq_along(trace_stages())) {
    message("  ", i, ". ", trace_stages()[[i]])
  }

  message("TRACE-Omics controlled perturbations")
  for (i in seq_along(trace_experiments())) {
    message("  ", i, ". ", trace_experiments()[[i]])
  }

  invisible(TRUE)
}

# ==============================================================================
# 4. Input table helpers
# ==============================================================================

trim_character <- function(x) {
  if (is.character(x)) trimws(x) else x
}

clean_colnames <- function(x) {
  colnames(x) <- trimws(as.character(colnames(x)))
  x
}

check_required_columns <- function(df, required, object_name = "data frame") {
  missing <- setdiff(required, names(df))

  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

read_count_table <- function(path, symbol_col = "symbol") {
  check_file_exists(path, "Count table")

  df <- safe_read_csv(path)
  df <- as.data.frame(df, check.names = FALSE)
  df <- clean_colnames(df)

  check_required_columns(df, symbol_col, "Count table")

  df[[symbol_col]] <- trimws(as.character(df[[symbol_col]]))
  df <- df[!is.na(df[[symbol_col]]) & nzchar(df[[symbol_col]]), , drop = FALSE]
  if (anyDuplicated(df[[symbol_col]])) {
    stop("Duplicated gene identifiers in count table: ", path, call. = FALSE)
  }

  df
}

read_metadata_table <- function(path,
                                sample_col = "sample_id",
                                cohort_col = "cohort",
                                outcome_col = "outcome") {
  check_file_exists(path, "Metadata file")

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    meta <- readRDS(path)
  } else {
    meta <- safe_read_csv(path)
  }

  meta <- as.data.frame(meta, check.names = FALSE)
  meta <- clean_colnames(meta)

  check_required_columns(
    meta,
    c(sample_col, cohort_col, outcome_col),
    "Metadata"
  )

  meta[[sample_col]] <- trimws(as.character(meta[[sample_col]]))
  meta[[cohort_col]] <- factor(trimws(as.character(meta[[cohort_col]])))
  meta[[outcome_col]] <- factor(trimws(as.character(meta[[outcome_col]])))

  if (anyDuplicated(meta[[sample_col]]) > 0) {
    stop("Metadata contains duplicated sample IDs.", call. = FALSE)
  }

  meta
}

merge_count_tables_by_symbol <- function(count_files, symbol_col = "symbol") {
  if (length(count_files) < 2) {
    stop("At least two count tables are required.", call. = FALSE)
  }

  count_list <- lapply(count_files, read_count_table, symbol_col = symbol_col)
  names(count_list) <- tools::file_path_sans_ext(basename(count_files))

  shared_genes <- Reduce(
    intersect,
    lapply(count_list, function(x) x[[symbol_col]])
  )

  if (length(shared_genes) == 0) {
    stop("No shared genes found across count tables.", call. = FALSE)
  }

  message("Shared genes across cohorts: ", length(shared_genes))

  count_list <- lapply(count_list, function(df) {
    df <- df[df[[symbol_col]] %in% shared_genes, , drop = FALSE]
    df <- df[order(df[[symbol_col]]), , drop = FALSE]
    df
  })

  merged <- Reduce(
    function(x, y) merge(x, y, by = symbol_col, all = FALSE, sort = FALSE),
    count_list
  )

  rownames(merged) <- make.unique(as.character(merged[[symbol_col]]))
  merged[[symbol_col]] <- NULL

  if (anyDuplicated(names(merged))) {
    stop("Duplicated sample identifiers across count tables.", call. = FALSE)
  }
  mat <- as.matrix(merged)
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (anyNA(mat)) {
    stop("Merged count matrix contains NA values.", call. = FALSE)
  }

  if (!all(is.finite(mat))) {
    stop("Merged count matrix contains non-finite values.", call. = FALSE)
  }

  if (any(mat < 0) || any(abs(mat - round(mat)) > 1e-8)) {
    stop("Count matrix must contain non-negative integer values; do not round normalised expression.", call. = FALSE)
  }
  mat <- round(mat)
  if (any(mat > .Machine$integer.max)) {
    stop("Count matrix has values exceeding R integer range.", call. = FALSE)
  }
  storage.mode(mat) <- "integer"

  mat
}

align_counts_metadata <- function(counts,
                                  metadata,
                                  sample_col = "sample_id",
                                  min_samples = 10) {
  counts <- as.matrix(counts)
  colnames(counts) <- trimws(as.character(colnames(counts)))

  metadata[[sample_col]] <- trimws(as.character(metadata[[sample_col]]))

  shared_samples <- intersect(colnames(counts), metadata[[sample_col]])

  if (length(shared_samples) < min_samples) {
    stop(
      "Too few matched samples between counts and metadata. ",
      "Matched samples: ", length(shared_samples),
      call. = FALSE
    )
  }

  message("Matched samples: ", length(shared_samples))
  message("Count samples dropped: ", length(setdiff(colnames(counts), shared_samples)))
  message("Metadata samples dropped: ", length(setdiff(metadata[[sample_col]], shared_samples)))

  counts <- counts[, shared_samples, drop = FALSE]
  metadata <- metadata[match(shared_samples, metadata[[sample_col]]), , drop = FALSE]

  if (!identical(colnames(counts), metadata[[sample_col]])) {
    stop("Sample alignment failed after matching.", call. = FALSE)
  }

  list(
    counts = counts,
    metadata = metadata
  )
}

# ==============================================================================
# 5. Matrix validation helpers
# ==============================================================================

coerce_numeric_matrix <- function(x, object_name = "matrix") {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  if (anyNA(x)) {
    stop(object_name, " contains NA values.", call. = FALSE)
  }

  if (!all(is.finite(x))) {
    stop(object_name, " contains non-finite values.", call. = FALSE)
  }

  x
}

remove_zero_variance_rows <- function(mat) {
  mat <- coerce_numeric_matrix(mat)

  vars <- apply(mat, 1, stats::var)
  keep <- is.finite(vars) & vars > 0

  mat[keep, , drop = FALSE]
}

align_matrix_list <- function(mat_list,
                              metadata,
                              sample_col = "sample_id",
                              min_samples = 10) {
  if (length(mat_list) == 0) {
    stop("mat_list is empty.", call. = FALSE)
  }

  mat_list <- lapply(mat_list, coerce_numeric_matrix)

  shared_samples <- Reduce(
    intersect,
    c(lapply(mat_list, colnames), list(metadata[[sample_col]]))
  )

  if (length(shared_samples) < min_samples) {
    stop(
      "Too few shared samples across matrices and metadata. ",
      "Shared samples: ", length(shared_samples),
      call. = FALSE
    )
  }

  mat_list <- lapply(mat_list, function(m) {
    m[, shared_samples, drop = FALSE]
  })

  metadata <- metadata[match(shared_samples, metadata[[sample_col]]), , drop = FALSE]

  for (nm in names(mat_list)) {
    if (!identical(colnames(mat_list[[nm]]), metadata[[sample_col]])) {
      stop("Sample alignment failed for matrix: ", nm, call. = FALSE)
    }
  }

  list(
    matrices = mat_list,
    metadata = metadata
  )
}

# ==============================================================================
# 6. Leakage-control and alignment audits
# ==============================================================================

assert_disjoint_samples <- function(development_ids,
                                    evaluation_ids,
                                    development_label = "development",
                                    evaluation_label = "evaluation") {
  development_ids <- unique(as.character(development_ids))
  evaluation_ids <- unique(as.character(evaluation_ids))

  overlap <- intersect(development_ids, evaluation_ids)

  if (length(overlap) > 0) {
    stop(
      "Sample leakage detected between ",
      development_label,
      " and ",
      evaluation_label,
      ": ",
      paste(head(overlap, 10), collapse = ", "),
      if (length(overlap) > 10) " ..." else "",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_feature_alignment <- function(reference_matrix,
                                     target_matrix,
                                     reference_label = "development",
                                     target_label = "evaluation") {
  reference_matrix <- as.matrix(reference_matrix)
  target_matrix <- as.matrix(target_matrix)

  reference_features <- rownames(reference_matrix)
  target_features <- rownames(target_matrix)

  if (is.null(reference_features) || is.null(target_features)) {
    stop("Both matrices require row names for feature-alignment auditing.", call. = FALSE)
  }

  if (!identical(reference_features, target_features)) {
    stop(
      "Feature order differs between ",
      reference_label,
      " and ",
      target_label,
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_binary_outcome <- function(metadata,
                                  outcome_col = "outcome",
                                  negative_class = "control",
                                  positive_class = "case") {
  check_required_columns(metadata, outcome_col, "Metadata")

  observed <- unique(as.character(metadata[[outcome_col]]))
  expected <- c(negative_class, positive_class)

  unexpected <- setdiff(observed, expected)

  if (length(unexpected) > 0) {
    stop(
      "Unexpected outcome level(s): ",
      paste(unexpected, collapse = ", "),
      call. = FALSE
    )
  }

  if (!all(expected %in% observed)) {
    stop(
      "Both outcome classes must be present: ",
      paste(expected, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_predictor_isolation <- function(feature_names,
                                       forbidden = c("sample_id", "participant_id", "study", "cohort", "cohort_group",
                                                  "outcome", "disease", "disease_group", "repository_accession", "partition")) {
  leaked <- intersect(as.character(feature_names), forbidden)

  if (length(leaked) > 0) {
    stop(
      "Non-molecular field(s) detected in predictor space: ",
      paste(leaked, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

write_sample_manifest <- function(sample_ids, path, label = "sample_id") {
  sample_ids <- unique(as.character(sample_ids))
  manifest <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  names(manifest)[[1]] <- label
  safe_write_csv(manifest, path)
  invisible(path)
}

# ==============================================================================
# 7. mlr3 task helpers
# ==============================================================================

is_mlr3_task <- function(x) {
  inherits(x, "Task")
}

read_mlr3_task <- function(path) {
  check_file_exists(path, "mlr3 task file")

  task <- readRDS(path)

  if (!inherits(task, "Task")) {
    stop("Object is not an mlr3 Task: ", path, call. = FALSE)
  }

  task
}

clean_classif_task <- function(task,
                               target = "outcome",
                               positive = "case",
                               remove_cols = c("cohort"),
                               clear_stratum = TRUE,
                               clear_group = TRUE) {
  if (!inherits(task, "TaskClassif")) {
    stop("Expected an mlr3 TaskClassif object.", call. = FALSE)
  }

  task <- task$clone(deep = TRUE)
  if (!positive %in% task$class_names) {
    stop("Requested positive class missing from task: ", positive, call. = FALSE)
  }
  task$positive <- positive

  removable <- intersect(remove_cols, task$feature_names)

  if (length(removable) > 0) {
    task$select(setdiff(task$feature_names, removable))
  }

  if (isTRUE(clear_stratum)) {
    task$col_roles$stratum <- character(0)
  }

  if (isTRUE(clear_group)) {
    task$col_roles$group <- character(0)
  }

  task
}

load_clean_classif_task <- function(path,
                                    positive = "case",
                                    remove_cols = c("cohort"),
                                    clear_stratum = TRUE,
                                    clear_group = TRUE) {
  task <- read_mlr3_task(path)

  clean_classif_task(
    task = task,
    positive = positive,
    remove_cols = remove_cols,
    clear_stratum = clear_stratum,
    clear_group = clear_group
  )
}

task_to_dataframe <- function(task) {
  if (!inherits(task, "Task")) {
    stop("Expected an mlr3 Task object.", call. = FALSE)
  }

  as.data.frame(task$data())
}

# ==============================================================================
# 8. Analytical label helpers
# ==============================================================================

# Ordered regex mapping implemented in base R to reduce CI-only dependencies.
# The first matching rule wins, preserving the legacy case_when() precedence.
trace_map_labels <- function(x, patterns, values, fallback = c("input", "NA", "None")) {
  fallback <- match.arg(fallback)
  x <- as.character(x)
  out <- if (fallback == "input") x else if (fallback == "NA")
    rep(NA_character_, length(x)) else rep("None", length(x))
  assigned <- rep(FALSE, length(x))
  for (i in seq_along(patterns)) {
    hit <- !assigned & !is.na(x) & grepl(patterns[[i]], x, perl = TRUE)
    out[hit] <- values[[i]]
    assigned[hit] <- TRUE
  }
  out
}

clean_representation_label <- function(x) {
  trace_map_labels(x,
    patterns = c(
      "VST_ComBatRef|ComBat[-_ ]?Ref",
      "ComBatSeq_log1p|ComBat[-_ ]?Seq.*log1p",
      "VST_fSVA|VST_SVA|fSVA|SVA",
      "VST_ComBat|VST [+] ComBat|ComBat",
      "VST_only|VST only|VST-only",
      "RAW_counts|Raw counts|raw counts",
      "log2[(]CPM[+]1[)]|logCPM|CPM"
    ),
    values = c("VST + ComBat-Ref", "ComBat-Seq + log1p", "VST + fSVA",
               "VST + ComBat", "VST only", "Raw counts", "log2(CPM+1)")
  )
}

clean_feature_selection_label <- function(x) {
  trace_map_labels(x,
    patterns = c("HybridA|Hybrid A", "HybridB|Hybrid B",
                 "information_gain|Info_Gain|Information gain",
                 "boruta|Boruta", "cmim|CMIM", "mrmr|mRMR",
                 "anova|ANOVA", "\\b(?:auc|AUC)\\b",
                 "importance|RF_Imp|RF importance", "variance|RF_Var|RF variance"),
    values = c("Hybrid A", "Hybrid B", "Information gain", "Boruta", "CMIM",
               "mRMR", "ANOVA", "AUC", "RF importance", "RF variance"),
    fallback = "NA"
  )
}

clean_class_handling_label <- function(x) {
  trace_map_labels(x,
    patterns = c("SMOTE", "Upsample|UpSample", "Downsample|DownSample",
                 "Weighted|Class weighted|class_weight",
                 "NBA|NoCBA|None|Unweighted"),
    values = c("SMOTE", "Upsample", "Downsample", "Class weighted", "None"),
    fallback = "None"
  )
}

clean_learner_label <- function(x) {
  trace_map_labels(x,
    patterns = c("Featureless", "Naive|naive_bayes", "Linear_SVM|Linear SVM",
                 "RBF_SVM|RBF SVM", "RandomForest|Random Forest|ranger",
                 "LightGBM|lightgbm", "TabPFN|tabpfn", "XGBoost|xgboost",
                 "C50|C5[.]0|C5"),
    values = c("Featureless", "Naive Bayes", "Linear SVM", "RBF SVM",
               "Random Forest", "LightGBM", "TabPFN", "XGBoost", "C5.0")
  )
}

make_configuration_label <- function(learner,
                                     representation,
                                     feature_selection,
                                     class_handling = "None") {
  paste(
    clean_learner_label(learner),
    clean_representation_label(representation),
    clean_feature_selection_label(feature_selection),
    clean_class_handling_label(class_handling),
    sep = " | "
  )
}

parse_task_id_trace <- function(task_id) {
  representation <- clean_representation_label(task_id)
  feature_selection <- clean_feature_selection_label(task_id)
  class_handling <- clean_class_handling_label(task_id)

  data.frame(
    task_id = task_id,
    representation = representation,
    feature_selection = feature_selection,
    class_handling = class_handling,
    HC = representation,
    FSF = feature_selection,
    CBA = class_handling,
    stringsAsFactors = FALSE
  )
}

# Backward-compatible aliases retained for earlier TRACE-Omics scripts.
clean_hc_label <- function(x) {
  clean_representation_label(x)
}

clean_fsf_label <- function(x) {
  clean_feature_selection_label(x)
}

clean_cba_label <- function(x) {
  clean_class_handling_label(x)
}

clean_mla_label <- function(x) {
  clean_learner_label(x)
}

# ==============================================================================
# 9. Prediction-table helpers
# ==============================================================================

find_prob_case_col <- function(dt) {
  nm <- names(dt)

  if ("prob.case" %in% nm) return("prob.case")
  if ("prob_case" %in% nm) return("prob_case")
  if ("prob.case." %in% nm) return("prob.case.")

  stop(
    "Could not find a case-probability column. Expected one of: ",
    "prob.case, prob_case, prob.case.",
    call. = FALSE
  )
}

standardise_prediction_table <- function(dt,
                                         positive_class = "case",
                                         negative_class = "control") {
  dt <- as.data.frame(dt, check.names = FALSE)

  required <- c("truth", "response")
  check_required_columns(dt, required, "Prediction table")

  prob_col <- find_prob_case_col(dt)

  dt$truth <- factor(as.character(dt$truth), levels = c(negative_class, positive_class))
  dt$response <- factor(as.character(dt$response), levels = c(negative_class, positive_class))
  dt$prob_case <- as.numeric(dt[[prob_col]])
  if (anyNA(dt$truth) || anyNA(dt$response) || anyNA(dt$prob_case) ||
      any(!is.finite(dt$prob_case)) || any(dt$prob_case < 0 | dt$prob_case > 1)) {
    stop("Prediction table has invalid labels or case probabilities.", call. = FALSE)
  }
  dt$y <- ifelse(as.character(dt$truth) == positive_class, 1, 0)

  dt
}

# ==============================================================================
# 10. Miscellaneous helpers
# ==============================================================================

message_section <- function(title) {
  line <- paste(rep("=", nchar(title) + 8), collapse = "")
  message("\n", line)
  message("=== ", title, " ===")
  message(line)
}

stop_if_empty <- function(x, object_name = "object") {
  if (length(x) == 0 || is.null(x)) {
    stop(object_name, " is empty.", call. = FALSE)
  }

  invisible(TRUE)
}

set_trace_seed <- function(seed = 7832) {
  set.seed(seed)
  invisible(seed)
}

make_trace_subseed <- function(seed = 7832, ...) {
  keys <- paste(..., collapse = "|")
  ints <- utf8ToInt(keys)

  if (length(ints) == 0) {
    return(as.integer(seed))
  }

  offset <- sum(ints * seq_along(ints))
  subseed <- as.integer((as.double(seed) + offset) %% .Machine$integer.max)

  if (subseed <= 0) {
    subseed <- as.integer(seed)
  }

  subseed
}

set_trace_subseed <- function(seed = 7832, ...) {
  subseed <- make_trace_subseed(seed, ...)
  set.seed(subseed)
  invisible(subseed)
}

# ==============================================================================
# End of TRACE-Omics R/utils_io.R
# =========================================