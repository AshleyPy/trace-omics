```r
# ==============================================================================
# TRACE-Omics
# Utility functions for input/output, task handling, configuration and validation
#
# File:
#   R/utils_io.R
#
# Purpose:
#   This file contains reusable helper functions used across the TRACE-Omics
#   workflow. The functions are intentionally lightweight and general so that the
#   main scripts remain easier to read, audit and adapt.
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
  df <- df[!duplicated(df[[symbol_col]]), , drop = FALSE]

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

  mat <- as.matrix(merged)
  storage.mode(mat) <- "numeric"

  if (anyNA(mat)) {
    stop("Merged count matrix contains NA values.", call. = FALSE)
  }

  if (!all(is.finite(mat))) {
    stop("Merged count matrix contains non-finite values.", call. = FALSE)
  }

  mat <- round(mat)
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
# 6. mlr3 task helpers
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
# 7. Label parsing helpers
# ==============================================================================

clean_hc_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "VST_ComBatRef|ComBatRef") ~ "VST + ComBat-Ref",
    stringr::str_detect(x, "VST_ComBat|ComBat") ~ "VST + ComBat",
    stringr::str_detect(x, "VST_SVA|SVA") ~ "VST + SVA",
    stringr::str_detect(x, "ComBatSeq_log1p") ~ "ComBat-Seq + log1p",
    stringr::str_detect(x, "VST_only|VST only") ~ "VST only",
    TRUE ~ as.character(x)
  )
}

clean_fsf_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "HybridA|Hybrid A") ~ "Hybrid A",
    stringr::str_detect(x, "HybridB|Hybrid B") ~ "Hybrid B",
    stringr::str_detect(x, "information_gain|Info_Gain") ~ "Info_Gain",
    stringr::str_detect(x, "boruta|Boruta") ~ "Boruta",
    stringr::str_detect(x, "cmim|CMIM") ~ "CMIM",
    stringr::str_detect(x, "mrmr|mRMR") ~ "mRMR",
    stringr::str_detect(x, "anova|ANOVA") ~ "ANOVA",
    stringr::str_detect(x, "auc|AUC") ~ "AUC",
    stringr::str_detect(x, "importance|RF_Imp") ~ "RF_Imp",
    stringr::str_detect(x, "variance|RF_Var") ~ "RF_Var",
    TRUE ~ NA_character_
  )
}

clean_cba_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "SMOTE$|SMOTE") ~ "SMOTE",
    stringr::str_detect(x, "Upsample$|Upsample") ~ "Upsample",
    stringr::str_detect(x, "Downsample$|Downsample") ~ "Downsample",
    stringr::str_detect(x, "NBA$|NBA") ~ "NBA",
    TRUE ~ "NBA"
  )
}

clean_mla_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "Featureless") ~ "Featureless",
    stringr::str_detect(x, "Naive") ~ "Naive Bayes",
    stringr::str_detect(x, "Linear_SVM|Linear SVM") ~ "Linear SVM",
    stringr::str_detect(x, "RBF_SVM|RBF SVM") ~ "RBF SVM",
    stringr::str_detect(x, "RandomForest|Random Forest") ~ "Random Forest",
    stringr::str_detect(x, "LightGBM") ~ "LightGBM",
    stringr::str_detect(x, "XGBoost") ~ "XGBoost",
    stringr::str_detect(x, "C50|C5") ~ "C50",
    TRUE ~ as.character(x)
  )
}

make_configuration_label <- function(MLA, HC, FSF, CBA = "NBA") {
  paste(
    clean_mla_label(MLA),
    clean_hc_label(HC),
    clean_fsf_label(FSF),
    clean_cba_label(CBA),
    sep = " | "
  )
}

parse_task_id_trace <- function(task_id) {
  data.frame(
    task_id = task_id,
    HC = clean_hc_label(task_id),
    FSF = clean_fsf_label(task_id),
    CBA = clean_cba_label(task_id),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 8. Prediction-table helpers
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
  dt$y <- ifelse(as.character(dt$truth) == positive_class, 1, 0)

  dt
}

# ==============================================================================
# 9. Miscellaneous helpers
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

# ==============================================================================
# End of R/utils_io.R
# ==============================================================================
```
