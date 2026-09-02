# ==============================================================================
# TRACE-Omics
# Stage 01: Representation specification
# ==============================================================================
#
# File:
#   scripts/01_representation_specification/01_representation_specification.R
#
# Framework:
#   Transcriptomic Robustness Assessment and Cross-cohort Evaluation for
#   Omics learning (TRACE-Omics)
#
# Purpose:
#   Reconstruct the complete representation-specification stage used by the
#   dissertation before feature selection and model fitting.
#
#   The script:
#     1. loads and audits the five disease-specific RNA-seq datasets;
#     2. harmonises cohort count tables on HGNC gene symbols;
#     3. creates or restores frozen development/evaluation manifests;
#     4. constructs development-only representation diagnostics;
#     5. prepares the representation inputs required by all four experiments:
#          Experiment 1 - sample-support perturbation;
#          Experiment 2 - disease-composition / LODO perturbation;
#          Experiment 3 - cohort-composition / LOCO perturbation;
#          Experiment 4 - matched LightGBM-TabPFN architecture perturbation;
#     6. saves fitted preprocessing states, matrices, manifests, diagnostics,
#        checksums and session information.
#
# IMPORTANT:
#   - All data-dependent transformations are fitted from the applicable
#     development/source observations.
#   - Evaluation observations enter only after the corresponding fitted state
#     has been estimated.
#   - ComBat-Ref and frozen ComBat-Seq projection are restricted to cohort
#     levels represented in the fitting population.
#   - Strict LOCO therefore uses VST/log2(CPM+1), and cross-disease LODO uses
#     VST + fSVA.
#
# Expected repository structure:
#
#   trace-omics/
#   |-- config/
#   |   `-- example_config.yml
#   |-- R/
#   |   `-- utils_io.R
#   |-- data/
#   |   `-- raw/
#   |       |-- TB/
#   |       |-- RA/
#   |       |-- UC/
#   |       |-- CD/
#   |       `-- T1D/
#   `-- scripts/
#       `-- 01_representation_specification/
#
# Example execution:
#
#   Rscript scripts/01_representation_specification/01_representation_specification.R
#
#   Rscript scripts/01_representation_specification/01_representation_specification.R \
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
    value <- if (length(pieces) >= 2) paste(pieces[-1], collapse = "=") else "true"

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
OVERWRITE_MANIFESTS <- parse_boolean(args$overwrite, default = FALSE)

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
# 1. LOAD TRACE-OMICS UTILITIES AND CONFIGURATION
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
    "readr",
    "ggplot2",
    "RColorBrewer",
    "cluster",
    "matrixStats"
  ),
  bioc = FALSE
)

require_packages(
  c(
    "DESeq2",
    "SummarizedExperiment",
    "sva",
    "edgeR",
    "neuroCombat"
  ),
  bioc = TRUE
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

config <- read_trace_config(CONFIG_PATH)

GLOBAL_SEED <- as.integer(get_config_value(config, "project.seed", 7832))
POSITIVE_CLASS <- as.character(get_config_value(config, "project.positive_class", "case"))
NEGATIVE_CLASS <- as.character(get_config_value(config, "project.negative_class", "control"))

SAMPLE_COL <- as.character(get_config_value(config, "analysis.sample_col", "sample_id"))
COHORT_COL <- as.character(get_config_value(config, "analysis.cohort_col", "cohort"))
OUTCOME_COL <- as.character(get_config_value(config, "analysis.outcome_col", "outcome"))

MIN_TOTAL_TRAINING_COUNT <- as.integer(get_config_value(config, "analysis.min_gene_count", 10))
DEFAULT_TOP_N <- as.integer(get_config_value(config, "analysis.top_n_features", 50))

MAX_SV <- 50L
MAX_PCA_K <- 50L
SILHOUETTE_NEAR_MAX <- as.numeric(
  get_config_value(
    config,
    "representation_diagnostics.silhouette.near_maximum_fraction",
    0.95
  )
)

EXPECTED_COMMON_CROSS_DISEASE_GENES <- as.integer(
  get_config_value(
    config,
    "disease_composition_experiment.common_feature_space.expected_gene_count",
    4450
  )
)

EXPECTED_MODEL_A_N <- as.integer(
  get_config_value(
    config,
    "disease_composition_experiment.model_a.expected_development_n",
    1413
  )
)

EXPECTED_MODEL_B_N <- as.integer(
  get_config_value(
    config,
    "disease_composition_experiment.model_b.expected_development_n",
    1736
  )
)

set_trace_seed(GLOBAL_SEED)


# ==============================================================================
# 2. PROJECT REGISTRY
# ==============================================================================

DISEASE_REGISTRY <- list(
  TB = list(
    expected_n = 465L,
    cohorts = c("GC6", "ACS", "Catalysis")
  ),
  RA = list(
    expected_n = 230L,
    cohorts = c("Sales", "Shchetynsky", "Guo")
  ),
  UC = list(
    expected_n = 1253L,
    cohorts = c("Wacker", "Halfvarson", "Mo")
  ),
  CD = list(
    expected_n = 495L,
    cohorts = c("Halfvarson", "Mo")
  ),
  T1D = list(
    expected_n = 55L,
    cohorts = c("Leal", "Speake")
  )
)

ALL_DISEASES <- names(DISEASE_REGISTRY)

if (!all(REQUESTED_DISEASES %in% ALL_DISEASES)) {
  stop(
    "Unknown disease identifier(s): ",
    paste(setdiff(REQUESTED_DISEASES, ALL_DISEASES), collapse = ", "),
    call. = FALSE
  )
}

RAW_ROOT <- dirname(get_config_value(config, "paths.raw_dir", "data/raw/TB"))
PROCESSED_ROOT <- dirname(get_config_value(config, "paths.processed_dir", "data/processed/TB"))

stage_example_path <- get_config_value(
  config,
  "paths.stage_results.representation_specification",
  "results/01_representation_specification/TB"
)

STAGE_ROOT <- dirname(stage_example_path)
AUDIT_ROOT <- file.path(STAGE_ROOT, "00_audits")
SESSION_ROOT <- file.path(STAGE_ROOT, "99_session_info")

make_dir(STAGE_ROOT)
make_dir(AUDIT_ROOT)
make_dir(SESSION_ROOT)

message_section("TRACE-Omics Stage 01: Representation specification")
message("Configuration: ", CONFIG_PATH)
message("Modes:         ", paste(RUN_MODES, collapse = ", "))
message("Diseases:      ", paste(REQUESTED_DISEASES, collapse = ", "))
message("Global seed:   ", GLOBAL_SEED)
message("Strict audits: ", STRICT_AUDITS)


# ==============================================================================
# 3. GENERAL REPRODUCIBILITY HELPERS
# ==============================================================================

audit_expected_value <- function(observed,
                                 expected,
                                 label,
                                 strict = STRICT_AUDITS) {
  if (is.na(expected) || is.null(expected)) {
    return(invisible(TRUE))
  }

  if (!identical(as.integer(observed), as.integer(expected))) {
    msg <- paste0(
      label,
      " audit differed from the dissertation expectation. Observed = ",
      observed,
      "; expected = ",
      expected,
      "."
    )

    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

write_md5_manifest <- function(paths, output_path) {
  paths <- unique(paths[file.exists(paths)])

  if (length(paths) == 0L) {
    return(invisible(NULL))
  }

  checksums <- tools::md5sum(paths)

  dt <- data.table(
    path = names(checksums),
    md5 = unname(checksums)
  )

  safe_write_csv(dt, output_path)

  invisible(dt)
}

save_config_snapshot <- function(config_path, output_dir) {
  make_dir(output_dir)

  target <- file.path(
    output_dir,
    paste0(
      "config_snapshot_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".yml"
    )
  )

  ok <- file.copy(config_path, target, overwrite = TRUE)

  if (!isTRUE(ok)) {
    warning("Configuration snapshot could not be written.", call. = FALSE)
  }

  invisible(target)
}

standardise_metadata <- function(metadata, disease_id) {
  metadata <- as.data.frame(metadata, check.names = FALSE)

  check_required_columns(
    metadata,
    c(SAMPLE_COL, COHORT_COL, OUTCOME_COL),
    paste0(disease_id, " metadata")
  )

  metadata[[SAMPLE_COL]] <- trimws(as.character(metadata[[SAMPLE_COL]]))
  metadata[[COHORT_COL]] <- factor(trimws(as.character(metadata[[COHORT_COL]])))
  metadata[[OUTCOME_COL]] <- factor(
    trimws(as.character(metadata[[OUTCOME_COL]])),
    levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)
  )

  if (anyDuplicated(metadata[[SAMPLE_COL]]) > 0L) {
    stop(disease_id, " metadata contains duplicated sample IDs.", call. = FALSE)
  }

  assert_binary_outcome(
    metadata,
    outcome_col = OUTCOME_COL,
    negative_class = NEGATIVE_CLASS,
    positive_class = POSITIVE_CLASS
  )

  metadata$sample_id <- metadata[[SAMPLE_COL]]
  metadata$cohort_group <- factor(as.character(metadata[[COHORT_COL]]))
  metadata$outcome <- factor(
    as.character(metadata[[OUTCOME_COL]]),
    levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)
  )
  metadata$disease <- disease_id

  metadata
}

load_disease_data <- function(disease_id) {
  raw_dir <- file.path(RAW_ROOT, disease_id)
  counts_dir <- file.path(raw_dir, "counts")

  metadata_candidates <- c(
    file.path(raw_dir, "metadata.rds"),
    file.path(raw_dir, "metadata.csv")
  )

  metadata_path <- metadata_candidates[file.exists(metadata_candidates)][1]

  if (is.na(metadata_path)) {
    stop(
      disease_id,
      " metadata file was not found. Expected metadata.rds or metadata.csv under ",
      raw_dir,
      call. = FALSE
    )
  }

  check_dir_exists(counts_dir, paste0(disease_id, " count directory"))

  count_files <- sort(
    list.files(
      counts_dir,
      pattern = "\\.csv$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  )

  if (length(count_files) < 2L) {
    stop(
      disease_id,
      " requires at least two cohort count tables in ",
      counts_dir,
      call. = FALSE
    )
  }

  message_section(paste0("Load and harmonise ", disease_id))

  counts <- merge_count_tables_by_symbol(
    count_files = count_files,
    symbol_col = "symbol"
  )

  metadata <- read_metadata_table(
    metadata_path,
    sample_col = SAMPLE_COL,
    cohort_col = COHORT_COL,
    outcome_col = OUTCOME_COL
  )

  metadata <- standardise_metadata(metadata, disease_id)

  aligned <- align_counts_metadata(
    counts = counts,
    metadata = metadata,
    sample_col = "sample_id",
    min_samples = 2L
  )

  counts <- aligned$counts
  metadata <- aligned$metadata

  counts <- counts[, metadata$sample_id, drop = FALSE]

  expected_n <- DISEASE_REGISTRY[[disease_id]]$expected_n

  audit_expected_value(
    observed = ncol(counts),
    expected = expected_n,
    label = paste0(disease_id, " final analytical population")
  )

  processed_dir <- file.path(PROCESSED_ROOT, disease_id, "representation_specification")
  manifest_dir <- file.path(PROCESSED_ROOT, disease_id, "manifests")

  make_dir(processed_dir)
  make_dir(manifest_dir)

  safe_save_rds(
    counts,
    file.path(processed_dir, paste0(disease_id, "_canonical_raw_counts.rds"))
  )

  safe_save_rds(
    metadata,
    file.path(processed_dir, paste0(disease_id, "_canonical_metadata.rds"))
  )

  write_md5_manifest(
    c(count_files, metadata_path),
    file.path(processed_dir, paste0(disease_id, "_input_md5.csv"))
  )

  audit <- data.table(
    disease = disease_id,
    genes = nrow(counts),
    samples = ncol(counts),
    controls = sum(metadata$outcome == NEGATIVE_CLASS),
    cases = sum(metadata$outcome == POSITIVE_CLASS),
    cohorts = length(unique(metadata$cohort_group))
  )

  safe_write_csv(
    audit,
    file.path(processed_dir, paste0(disease_id, "_canonical_input_audit.csv"))
  )

  list(
    disease_id = disease_id,
    counts = counts,
    metadata = metadata,
    raw_dir = raw_dir,
    count_files = count_files,
    metadata_path = metadata_path,
    processed_dir = processed_dir,
    manifest_dir = manifest_dir
  )
}


# ==============================================================================
# 4. STRATIFIED ALLOCATION AND MANIFEST HELPERS
# ==============================================================================

allocate_integer_target <- function(counts_by_stratum,
                                    target,
                                    retain_at_least_one = FALSE,
                                    give_each_if_possible = TRUE) {
  stratum_names <- names(counts_by_stratum)

  counts_by_stratum <- as.integer(counts_by_stratum)
  names(counts_by_stratum) <- stratum_names

  if (target < 0L) {
    stop("Allocation target must be non-negative.", call. = FALSE)
  }

  if (sum(counts_by_stratum) < target) {
    stop("Allocation target exceeds available observations.", call. = FALSE)
  }

  eligible <- counts_by_stratum > 0L

  allocation <- integer(length(counts_by_stratum))
  names(allocation) <- names(counts_by_stratum)

  capacity <- counts_by_stratum

  if (retain_at_least_one) {
    capacity <- pmax(
      counts_by_stratum - as.integer(eligible),
      0L
    )
  }

  if (sum(capacity) < target) {
    stop(
      "Requested allocation cannot retain at least one observation in every non-empty stratum.",
      call. = FALSE
    )
  }

  if (
    isTRUE(give_each_if_possible) &&
      target >= sum(eligible & capacity > 0L)
  ) {
    initial <- which(eligible & capacity > 0L)
    allocation[initial] <- 1L
  }

  allocation <- pmin(allocation, capacity)
  remaining <- target - sum(allocation)

  while (remaining > 0L) {
    available <- which(allocation < capacity)

    if (length(available) == 0L) {
      stop("Allocation ended before reaching the requested target.", call. = FALSE)
    }

    residual_capacity <- capacity[available] - allocation[available]
    weights <- counts_by_stratum[available]

    if (sum(weights) <= 0L) {
      weights <- residual_capacity
    }

    raw_add <- remaining * weights / sum(weights)
    floor_add <- pmin(floor(raw_add), residual_capacity)

    if (sum(floor_add) > 0L) {
      allocation[available] <- allocation[available] + floor_add
      remaining <- target - sum(allocation)

      if (remaining <= 0L) {
        break
      }
    }

    available <- which(allocation < capacity)

    if (length(available) == 0L) {
      break
    }

    fractional <- vapply(
      available,
      function(i) {
        raw <- target * counts_by_stratum[[i]] / sum(counts_by_stratum)
        raw - floor(raw)
      },
      numeric(1)
    )

    order_add <- order(
      fractional,
      counts_by_stratum[available],
      decreasing = TRUE
    )

    for (ordered_index in order_add) {
      if (remaining <= 0L) {
        break
      }

      stratum_index <- available[[ordered_index]]

      if (allocation[[stratum_index]] < capacity[[stratum_index]]) {
        allocation[[stratum_index]] <- allocation[[stratum_index]] + 1L
        remaining <- remaining - 1L
      }
    }
  }

  if (sum(allocation) != target) {
    stop("Final allocation differs from the requested target.", call. = FALSE)
  }

  allocation
}

select_ids_by_outcome_and_cohort <- function(metadata,
                                             outcome_value,
                                             target,
                                             seed,
                                             retain_at_least_one = FALSE) {
  current <- metadata[metadata$outcome == outcome_value, , drop = FALSE]

  cohort_counts <- table(current$cohort_group)

  allocation <- allocate_integer_target(
    counts_by_stratum = cohort_counts,
    target = target,
    retain_at_least_one = retain_at_least_one,
    give_each_if_possible = TRUE
  )

  selected_ids <- character(0)
  allocation_rows <- list()

  for (cohort_index in seq_along(allocation)) {
    cohort_name <- names(allocation)[[cohort_index]]
    requested <- allocation[[cohort_name]]

    cohort_ids <- current$sample_id[
      as.character(current$cohort_group) == cohort_name
    ]

    set_trace_subseed(
      seed,
      "outcome",
      outcome_value,
      "cohort",
      cohort_name,
      "selection"
    )

    chosen <- if (requested > 0L) {
      sample(
        cohort_ids,
        size = requested,
        replace = FALSE
      )
    } else {
      character(0)
    }

    selected_ids <- c(selected_ids, chosen)

    allocation_rows[[cohort_index]] <- data.table(
      outcome = outcome_value,
      cohort = cohort_name,
      available = length(cohort_ids),
      selected = requested,
      retained = length(cohort_ids) - requested
    )
  }

  list(
    ids = selected_ids,
    allocation = rbindlist(allocation_rows, fill = TRUE)
  )
}

create_general_development_manifest <- function(metadata,
                                                disease_id,
                                                manifest_path,
                                                development_fraction = 0.80,
                                                seed = GLOBAL_SEED,
                                                overwrite = FALSE) {
  if (file.exists(manifest_path) && !isTRUE(overwrite)) {
    manifest <- safe_read_csv(manifest_path)
    manifest <- as.data.frame(manifest, check.names = FALSE)

    check_required_columns(
      manifest,
      c("sample_id", "role"),
      paste0(disease_id, " partition manifest")
    )

    if (!setequal(manifest$sample_id, metadata$sample_id)) {
      stop(
        disease_id,
        " saved partition manifest contains a different sample population.",
        call. = FALSE
      )
    }

    return(manifest)
  }

  holdout_fraction <- 1 - development_fraction
  holdout_target <- as.integer(round(nrow(metadata) * holdout_fraction))

  strata <- interaction(
    metadata$cohort_group,
    metadata$outcome,
    drop = TRUE,
    lex.order = TRUE
  )

  counts_by_stratum <- table(strata)

  allocation <- allocate_integer_target(
    counts_by_stratum = counts_by_stratum,
    target = holdout_target,
    retain_at_least_one = TRUE,
    give_each_if_possible = TRUE
  )

  holdout_ids <- character(0)

  for (stratum_index in seq_along(allocation)) {
    stratum_name <- names(allocation)[[stratum_index]]
    requested <- allocation[[stratum_name]]

    ids <- metadata$sample_id[as.character(strata) == stratum_name]

    set_trace_subseed(
      seed,
      disease_id,
      "general_split",
      stratum_name
    )

    selected <- if (requested > 0L) {
      sample(ids, requested, replace = FALSE)
    } else {
      character(0)
    }

    holdout_ids <- c(holdout_ids, selected)
  }

  role <- ifelse(
    metadata$sample_id %in% holdout_ids,
    "held_out",
    "development"
  )

  manifest <- data.frame(
    disease = disease_id,
    sample_id = metadata$sample_id,
    cohort = as.character(metadata$cohort_group),
    outcome = as.character(metadata$outcome),
    role = role,
    stringsAsFactors = FALSE
  )

  assert_disjoint_samples(
    development_ids = manifest$sample_id[manifest$role == "development"],
    evaluation_ids = manifest$sample_id[manifest$role == "held_out"]
  )

  safe_write_csv(manifest, manifest_path)

  manifest
}

subset_from_manifest <- function(disease_object, manifest) {
  development_ids <- manifest$sample_id[manifest$role == "development"]
  held_out_ids <- manifest$sample_id[manifest$role == "held_out"]

  metadata <- disease_object$metadata
  counts <- disease_object$counts

  development_metadata <- metadata[
    match(development_ids, metadata$sample_id),
    ,
    drop = FALSE
  ]

  held_out_metadata <- metadata[
    match(held_out_ids, metadata$sample_id),
    ,
    drop = FALSE
  ]

  development_counts <- counts[
    ,
    development_metadata$sample_id,
    drop = FALSE
  ]

  held_out_counts <- counts[
    ,
    held_out_metadata$sample_id,
    drop = FALSE
  ]

  assert_disjoint_samples(
    development_metadata$sample_id,
    held_out_metadata$sample_id
  )

  list(
    development_counts = development_counts,
    development_metadata = development_metadata,
    held_out_counts = held_out_counts,
    held_out_metadata = held_out_metadata
  )
}

make_stratified_fold_map <- function(metadata,
                                     folds = 5L,
                                     seed = GLOBAL_SEED) {
  if (folds < 2L) {
    stop("At least two folds are required.", call. = FALSE)
  }

  strata <- interaction(
    metadata$cohort_group,
    metadata$outcome,
    drop = TRUE,
    lex.order = TRUE
  )

  fold_id <- integer(nrow(metadata))

  for (stratum_name in levels(strata)) {
    rows <- which(strata == stratum_name)

    set_trace_subseed(
      seed,
      "fold_map",
      stratum_name,
      folds
    )

    shuffled <- sample(rows, length(rows), replace = FALSE)
    assigned <- rep(seq_len(folds), length.out = length(rows))

    fold_id[shuffled] <- assigned
  }

  data.frame(
    sample_id = metadata$sample_id,
    disease = metadata$disease,
    cohort = as.character(metadata$cohort_group),
    outcome = as.character(metadata$outcome),
    fold = fold_id,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 5. FROZEN VST
# ==============================================================================

calculate_poscounts_geo_means <- function(count_matrix) {
  positive_number <- rowSums(count_matrix > 0)

  if (any(positive_number == 0L)) {
    stop(
      "At least one retained gene is zero in every fitting observation.",
      call. = FALSE
    )
  }

  log_counts <- log(count_matrix)
  log_counts[!is.finite(log_counts)] <- 0

  geo_means <- exp(
    rowSums(log_counts) /
      positive_number
  )

  names(geo_means) <- rownames(count_matrix)

  if (
    anyNA(geo_means) ||
      !all(is.finite(geo_means)) ||
      any(geo_means <= 0)
  ) {
    stop("Invalid fitting-set geometric means.", call. = FALSE)
  }

  geo_means
}

make_training_design <- function(metadata,
                                 include_cohort = TRUE,
                                 include_outcome = TRUE) {
  terms <- character(0)

  if (isTRUE(include_cohort) && length(unique(metadata$cohort_group)) > 1L) {
    terms <- c(terms, "cohort_group")
  }

  if (isTRUE(include_outcome) && length(unique(metadata$outcome)) > 1L) {
    terms <- c(terms, "outcome")
  }

  if (length(terms) == 0L) {
    return(~ 1)
  }

  stats::as.formula(
    paste("~", paste(terms, collapse = " + "))
  )
}

fit_apply_frozen_vst <- function(train_counts,
                                 new_counts = NULL,
                                 train_metadata,
                                 min_total_count = MIN_TOTAL_TRAINING_COUNT,
                                 min_features = DEFAULT_TOP_N,
                                 include_cohort_in_design = TRUE) {
  train_counts <- train_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]

  training_total_counts <- rowSums(train_counts)

  keep_gene <- (
    is.finite(training_total_counts) &
      training_total_counts >= min_total_count
  )

  if (sum(keep_gene) < min_features) {
    stop(
      "Fewer than ",
      min_features,
      " genes survived the fitting-only count filter.",
      call. = FALSE
    )
  }

  train_counts_filtered <- train_counts[
    keep_gene,
    ,
    drop = FALSE
  ]

  new_counts_filtered <- if (!is.null(new_counts)) {
    new_counts[
      rownames(train_counts_filtered),
      ,
      drop = FALSE
    ]
  } else {
    NULL
  }

  train_coldata <- data.frame(
    cohort_group = factor(as.character(train_metadata$cohort_group)),
    outcome = factor(
      as.character(train_metadata$outcome),
      levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)
    ),
    row.names = train_metadata$sample_id,
    check.names = FALSE
  )

  design_formula <- make_training_design(
    train_coldata,
    include_cohort = include_cohort_in_design,
    include_outcome = TRUE
  )

  training_geo_means <- calculate_poscounts_geo_means(
    train_counts_filtered
  )

  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData = train_counts_filtered,
    colData = train_coldata,
    design = design_formula
  )

  dds_train <- DESeq2::estimateSizeFactors(
    dds_train,
    type = "poscounts",
    geoMeans = training_geo_means
  )

  dds_train <- tryCatch(
    DESeq2::estimateDispersions(
      dds_train,
      fitType = "parametric",
      quiet = TRUE
    ),
    error = function(e) {
      DESeq2::estimateDispersions(
        dds_train,
        fitType = "local",
        quiet = TRUE
      )
    }
  )

  training_dispersion_function <- DESeq2::dispersionFunction(
    dds_train
  )

  vst_train_object <- DESeq2::varianceStabilizingTransformation(
    dds_train,
    blind = FALSE
  )

  train_vst <- as.matrix(
    SummarizedExperiment::assay(
      vst_train_object
    )
  )

  new_vst <- NULL

  if (!is.null(new_counts_filtered)) {
    new_ids <- colnames(new_counts_filtered)

    new_coldata <- data.frame(
      intercept = rep(1L, length(new_ids)),
      row.names = new_ids
    )

    dds_new <- DESeq2::DESeqDataSetFromMatrix(
      countData = new_counts_filtered,
      colData = new_coldata,
      design = ~ 1
    )

    dds_new <- DESeq2::estimateSizeFactors(
      dds_new,
      type = "poscounts",
      geoMeans = training_geo_means
    )

    DESeq2::dispersionFunction(dds_new) <- training_dispersion_function

    vst_new_object <- DESeq2::varianceStabilizingTransformation(
      dds_new,
      blind = FALSE
    )

    new_vst <- as.matrix(
      SummarizedExperiment::assay(
        vst_new_object
      )
    )

    assert_feature_alignment(
      train_vst,
      new_vst,
      reference_label = "VST fitting matrix",
      target_label = "VST evaluation matrix"
    )
  }

  list(
    train = train_vst,
    new = new_vst,
    state = list(
      representation = "VST_only",
      genes = rownames(train_counts_filtered),
      geo_means = training_geo_means,
      dispersion_function = training_dispersion_function,
      training_design = deparse(design_formula),
      min_total_count = min_total_count
    )
  )
}


# ==============================================================================
# 6. FROZEN VST + fSVA
# ==============================================================================

fit_apply_vst_fsva <- function(train_counts,
                               new_counts = NULL,
                               train_metadata,
                               seed,
                               include_cohort_in_vst_design = TRUE) {
  vst <- fit_apply_frozen_vst(
    train_counts = train_counts,
    new_counts = new_counts,
    train_metadata = train_metadata,
    include_cohort_in_design = include_cohort_in_vst_design
  )

  train_vst <- vst$train
  new_vst <- vst$new

  training_variance <- matrixStats::rowVars(train_vst)

  sva_gene_keep <- (
    is.finite(training_variance) &
      training_variance > sqrt(.Machine$double.eps)
  )

  train_sva_input <- train_vst[
    sva_gene_keep,
    ,
    drop = FALSE
  ]

  new_sva_input <- if (!is.null(new_vst)) {
    new_vst[
      sva_gene_keep,
      ,
      drop = FALSE
    ]
  } else {
    train_sva_input
  }

  if (nrow(train_sva_input) < DEFAULT_TOP_N) {
    stop(
      "Fewer than ",
      DEFAULT_TOP_N,
      " variable genes were available for fSVA.",
      call. = FALSE
    )
  }

  train_coldata <- data.frame(
    outcome = factor(
      as.character(train_metadata$outcome),
      levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)
    ),
    row.names = train_metadata$sample_id
  )

  mod_full <- stats::model.matrix(
    ~ outcome,
    data = train_coldata
  )

  mod_null <- stats::model.matrix(
    ~ 1,
    data = train_coldata
  )

  set.seed(seed)

  n_sv_be <- tryCatch(
    as.integer(
      sva::num.sv(
        dat = train_sva_input,
        mod = mod_full,
        method = "be"
      )
    ),
    error = function(e) NA_integer_
  )

  n_sv_leek <- tryCatch(
    as.integer(
      sva::num.sv(
        dat = train_sva_input,
        mod = mod_full,
        method = "leek"
      )
    ),
    error = function(e) NA_integer_
  )

  if (!is.na(n_sv_be) && is.finite(n_sv_be)) {
    n_sv_requested <- n_sv_be
  } else if (!is.na(n_sv_leek) && is.finite(n_sv_leek)) {
    n_sv_requested <- n_sv_leek
  } else {
    n_sv_requested <- 0L
  }

  max_sv_by_model <- max(
    0L,
    ncol(train_sva_input) - qr(mod_full)$rank - 1L
  )

  max_sv_by_sample_cap <- floor(
    ncol(train_sva_input) / 3L
  )

  max_allowed_sv <- max(
    0L,
    min(
      MAX_SV,
      max_sv_by_model,
      max_sv_by_sample_cap,
      nrow(train_sva_input) - 1L
    )
  )

  n_sv_used <- max(
    0L,
    min(
      n_sv_requested,
      max_allowed_sv
    )
  )

  train_fsva <- train_vst
  new_fsva <- new_vst
  sva_fit <- NULL

  if (n_sv_used > 0L) {
    set_trace_subseed(
      seed,
      "fsva_fit",
      n_sv_used
    )

    sva_fit <- sva::sva(
      dat = train_sva_input,
      mod = mod_full,
      mod0 = mod_null,
      n.sv = n_sv_used
    )

    fsva_result <- sva::fsva(
      dbdat = train_sva_input,
      mod = mod_full,
      sv = sva_fit,
      newdat = new_sva_input,
      method = "fast"
    )

    train_corrected_subset <- as.matrix(fsva_result$db)
    dimnames(train_corrected_subset) <- dimnames(train_sva_input)

    train_fsva[sva_gene_keep, ] <- train_corrected_subset

    if (!is.null(new_vst)) {
      new_corrected_subset <- as.matrix(fsva_result$new)
      dimnames(new_corrected_subset) <- dimnames(new_sva_input)
      new_fsva[sva_gene_keep, ] <- new_corrected_subset
    }
  }

  if (anyNA(train_fsva) || !all(is.finite(train_fsva))) {
    stop("Invalid VST + fSVA fitting matrix.", call. = FALSE)
  }

  if (
    !is.null(new_fsva) &&
      (anyNA(new_fsva) || !all(is.finite(new_fsva)))
  ) {
    stop("Invalid VST + fSVA evaluation matrix.", call. = FALSE)
  }

  list(
    train = train_fsva,
    new = new_fsva,
    state = c(
      vst$state,
      list(
        representation = "VST_fSVA",
        sva_genes = rownames(train_vst)[sva_gene_keep],
        nSV_BE = n_sv_be,
        nSV_Leek = n_sv_leek,
        nSV_used = n_sv_used,
        sva_object = sva_fit
      )
    )
  )
}


# ==============================================================================
# 7. FROZEN VST + COMBAT / COMBAT-REF
# ==============================================================================

fit_apply_vst_combat <- function(train_counts,
                                 new_counts = NULL,
                                 train_metadata,
                                 new_metadata = NULL,
                                 reference_anchored = FALSE) {
  vst <- fit_apply_frozen_vst(
    train_counts = train_counts,
    new_counts = new_counts,
    train_metadata = train_metadata,
    include_cohort_in_design = TRUE
  )

  expr_train <- vst$train
  expr_new <- vst$new

  train_batches <- as.character(train_metadata$cohort_group)

  if (length(unique(train_batches)) < 2L) {
    stop("ComBat requires at least two fitting cohorts.", call. = FALSE)
  }

  if (any(table(train_batches) < 2L)) {
    stop(
      "ComBat requires at least two fitting observations per cohort.",
      call. = FALSE
    )
  }

  if (!is.null(expr_new)) {
    if (is.null(new_metadata)) {
      stop(
        "new_metadata is required when projecting ComBat to evaluation observations.",
        call. = FALSE
      )
    }

    new_batches <- as.character(new_metadata$cohort_group)

    unseen_batches <- setdiff(
      unique(new_batches),
      unique(train_batches)
    )

    if (length(unseen_batches) > 0L) {
      stop(
        "Frozen ComBat cannot transform cohort(s) absent from fitting data: ",
        paste(unseen_batches, collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    new_batches <- NULL
  }

  within_batch_residual_ss <- numeric(nrow(expr_train))

  for (batch_level in unique(train_batches)) {
    batch_columns <- which(train_batches == batch_level)

    batch_matrix <- expr_train[
      ,
      batch_columns,
      drop = FALSE
    ]

    batch_centred <- sweep(
      batch_matrix,
      MARGIN = 1L,
      STATS = rowMeans(batch_matrix),
      FUN = "-"
    )

    within_batch_residual_ss <- (
      within_batch_residual_ss +
        rowSums(batch_centred^2)
    )
  }

  combat_gene_keep <- (
    is.finite(within_batch_residual_ss) &
      within_batch_residual_ss > sqrt(.Machine$double.eps)
  )

  if (!any(combat_gene_keep)) {
    stop(
      "No genes contained sufficient within-cohort fitting variation for ComBat.",
      call. = FALSE
    )
  }

  expr_train_input <- expr_train[
    combat_gene_keep,
    ,
    drop = FALSE
  ]

  ref_level <- if (isTRUE(reference_anchored)) {
    names(sort(table(train_batches), decreasing = TRUE))[[1]]
  } else {
    NULL
  }

  fit_train <- neuroCombat::neuroCombat(
    dat = expr_train_input,
    batch = train_batches,
    mod = NULL,
    eb = TRUE,
    parametric = TRUE,
    mean.only = FALSE,
    ref.batch = ref_level,
    verbose = FALSE
  )

  adjusted_train <- expr_train

  adjusted_train[
    combat_gene_keep,

  ] <- as.matrix(fit_train$dat.combat)

  adjusted_new <- NULL

  if (!is.null(expr_new)) {
    expr_new_input <- expr_new[
      combat_gene_keep,
      ,
      drop = FALSE
    ]

    fit_new <- neuroCombat::neuroCombatFromTraining(
      dat = expr_new_input,
      batch = new_batches,
      estimates = fit_train$estimates,
      mod = NULL,
      verbose = FALSE
    )

    adjusted_new <- expr_new

    adjusted_new[
      combat_gene_keep,

    ] <- as.matrix(fit_new$dat.combat)

    dimnames(adjusted_new) <- dimnames(expr_new)
  }

  representation_id <- if (isTRUE(reference_anchored)) {
    "VST_ComBatRef"
  } else {
    "VST_ComBat"
  }

  list(
    train = adjusted_train,
    new = adjusted_new,
    state = c(
      vst$state,
      list(
        representation = representation_id,
        reference_cohort = ref_level,
        combat_genes = rownames(expr_train)[combat_gene_keep],
        combat_estimates = fit_train$estimates
      )
    )
  )
}


# ==============================================================================
# 8. STANDARD COMBAT-SEQ + LOG1P FOR DEVELOPMENT DIAGNOSTICS
# ==============================================================================

fit_standard_combatseq_diagnostic <- function(train_counts,
                                              train_metadata) {
  training_total_counts <- rowSums(train_counts)

  keep_gene <- (
    is.finite(training_total_counts) &
      training_total_counts >= MIN_TOTAL_TRAINING_COUNT
  )

  counts_filtered <- train_counts[
    keep_gene,
    ,
    drop = FALSE
  ]

  adjusted_counts <- sva::ComBat_seq(
    counts = counts_filtered,
    batch = as.character(train_metadata$cohort_group),
    group = as.character(train_metadata$outcome),
    full_mod = TRUE,
    shrink = FALSE,
    shrink.disp = FALSE
  )

  transformed <- log1p(adjusted_counts)

  if (anyNA(transformed) || !all(is.finite(transformed))) {
    stop(
      "Invalid development-only ComBat-Seq + log1p representation.",
      call. = FALSE
    )
  }

  list(
    train = transformed,
    new = NULL,
    state = list(
      representation = "ComBatSeq_log1p",
      implementation = "sva::ComBat_seq development diagnostic",
      genes = rownames(counts_filtered),
      batch_variable = "cohort_group",
      biological_group = "outcome",
      full_mod = TRUE,
      shrink = FALSE,
      shrink_disp = FALSE
    )
  )
}


# ==============================================================================
# 9. FROZEN TRAINING-ESTIMATED COMBAT-SEQ + LOG1P
#    Used for Experiment 1 projection when evaluation cohort levels occur in
#    the fitting population.
# ==============================================================================

vec_to_matrix <- function(vec, n_columns) {
  matrix(
    rep(vec, n_columns),
    ncol = n_columns,
    byrow = FALSE
  )
}

combatseq_match_quantiles <- function(counts_sub,
                                      old_mu,
                                      old_phi,
                                      new_mu,
                                      new_phi) {
  adjusted <- matrix(
    NA_real_,
    nrow = nrow(counts_sub),
    ncol = ncol(counts_sub),
    dimnames = dimnames(counts_sub)
  )

  for (gene_index in seq_len(nrow(counts_sub))) {
    for (sample_index in seq_len(ncol(counts_sub))) {
      observed <- counts_sub[
        gene_index,
        sample_index
      ]

      if (observed <= 1L) {
        adjusted[
          gene_index,
          sample_index
        ] <- observed

        next
      }

      probability <- stats::pnbinom(
        observed - 1L,
        mu = old_mu[
          gene_index,
          sample_index
        ],
        size = 1 / old_phi[[gene_index]]
      )

      if (
        !is.finite(probability) ||
          abs(probability - 1) < 1e-4
      ) {
        adjusted[
          gene_index,
          sample_index
        ] <- observed
      } else {
        adjusted[
          gene_index,
          sample_index
        ] <- 1 +
          stats::qnbinom(
            probability,
            mu = new_mu[
              gene_index,
              sample_index
            ],
            size = 1 / new_phi[[gene_index]]
          )
      }
    }
  }

  adjusted
}

fit_frozen_combatseq <- function(train_counts,
                                 batch_train) {
  training_total_counts <- rowSums(train_counts)

  keep_initial <- (
    is.finite(training_total_counts) &
      training_total_counts >= MIN_TOTAL_TRAINING_COUNT
  )

  train_counts <- train_counts[
    keep_initial,
    ,
    drop = FALSE
  ]

  batch_train <- droplevels(
    factor(as.character(batch_train))
  )

  if (any(table(batch_train) <= 1L)) {
    stop(
      "Frozen ComBat-Seq requires at least two fitting observations in every cohort.",
      call. = FALSE
    )
  }

  keep_by_batch <- lapply(
    levels(batch_train),
    function(batch_level) {
      which(
        rowSums(
          train_counts[
            ,
            batch_train == batch_level,
            drop = FALSE
          ]
        ) > 0L
      )
    }
  )

  gene_keep_index <- Reduce(
    intersect,
    keep_by_batch
  )

  if (length(gene_keep_index) < DEFAULT_TOP_N) {
    stop(
      "Fewer than ",
      DEFAULT_TOP_N,
      " genes had non-zero counts in every fitting cohort for frozen ComBat-Seq.",
      call. = FALSE
    )
  }

  gene_remove_index <- setdiff(
    seq_len(nrow(train_counts)),
    gene_keep_index
  )

  counts_fit <- train_counts[
    gene_keep_index,
    ,
    drop = FALSE
  ]

  dge_object <- edgeR::DGEList(
    counts = counts_fit
  )

  batch_levels <- levels(batch_train)
  n_batches <- length(batch_levels)

  batch_indices <- lapply(
    batch_levels,
    function(batch_level) {
      which(batch_train == batch_level)
    }
  )

  batch_sizes <- vapply(
    batch_indices,
    length,
    integer(1)
  )

  n_samples <- sum(batch_sizes)

  batch_design <- stats::model.matrix(
    ~ -1 + batch_train
  )

  colnames(batch_design) <- batch_levels

  intercept_design <- matrix(
    1,
    nrow = n_samples,
    ncol = 1L
  )

  common_dispersion <- vapply(
    seq_len(n_batches),
    function(batch_index) {
      current_columns <- batch_indices[[batch_index]]

      edgeR::estimateGLMCommonDisp(
        counts_fit[
          ,
          current_columns,
          drop = FALSE
        ],
        design = NULL,
        subset = nrow(counts_fit)
      )
    },
    numeric(1)
  )

  gene_dispersion <- lapply(
    seq_len(n_batches),
    function(batch_index) {
      current_columns <- batch_indices[[batch_index]]

      edgeR::estimateGLMTagwiseDisp(
        counts_fit[
          ,
          current_columns,
          drop = FALSE
        ],
        design = intercept_design[
          current_columns,
          ,
          drop = FALSE
        ],
        dispersion = common_dispersion[[batch_index]],
        prior.df = 0
      )
    }
  )

  phi_matrix <- matrix(
    NA_real_,
    nrow = nrow(counts_fit),
    ncol = ncol(counts_fit)
  )

  for (batch_index in seq_len(n_batches)) {
    phi_matrix[
      ,
      batch_indices[[batch_index]]
    ] <- vec_to_matrix(
      gene_dispersion[[batch_index]],
      batch_sizes[[batch_index]]
    )
  }

  first_glm <- edgeR::glmFit(
    dge_object,
    design = batch_design,
    dispersion = phi_matrix,
    prior.count = 1e-4
  )

  alpha_g <- as.numeric(
    first_glm$coefficients[
      ,
      seq_len(n_batches),
      drop = FALSE
    ] %*%
      matrix(
        batch_sizes / n_samples,
        ncol = 1L
      )
  )

  names(alpha_g) <- rownames(counts_fit)

  training_offset <- edgeR::getOffset(
    dge_object
  )

  new_offset <- t(
    vec_to_matrix(
      training_offset,
      nrow(counts_fit)
    )
  ) +
    vec_to_matrix(
      alpha_g,
      ncol(counts_fit)
    )

  second_glm <- edgeR::glmFit(
    counts_fit,
    design = batch_design,
    dispersion = phi_matrix,
    offset = new_offset,
    prior.count = 1e-4
  )

  gamma_star <- second_glm$coefficients[
    ,
    seq_len(n_batches),
    drop = FALSE
  ]

  colnames(gamma_star) <- batch_levels
  rownames(gamma_star) <- rownames(counts_fit)

  phi_hat <- do.call(
    cbind,
    gene_dispersion
  )

  colnames(phi_hat) <- batch_levels
  rownames(phi_hat) <- rownames(counts_fit)

  phi_star <- rowMeans(phi_hat)

  mu_hat_train <- second_glm$fitted.values

  adjusted_fit_counts <- matrix(
    NA_real_,
    nrow = nrow(counts_fit),
    ncol = ncol(counts_fit),
    dimnames = dimnames(counts_fit)
  )

  for (batch_index in seq_len(n_batches)) {
    current_columns <- batch_indices[[batch_index]]

    new_mu <- exp(
      log(
        mu_hat_train[
          ,
          current_columns,
          drop = FALSE
        ]
      ) -
        vec_to_matrix(
          gamma_star[
            ,
            batch_index
          ],
          length(current_columns)
        )
    )

    adjusted_fit_counts[
      ,
      current_columns
    ] <- combatseq_match_quantiles(
      counts_sub = counts_fit[
        ,
        current_columns,
        drop = FALSE
      ],
      old_mu = mu_hat_train[
        ,
        current_columns,
        drop = FALSE
      ],
      old_phi = phi_hat[
        ,
        batch_index
      ],
      new_mu = new_mu,
      new_phi = phi_star
    )
  }

  adjusted_train <- matrix(
    NA_real_,
    nrow = nrow(train_counts),
    ncol = ncol(train_counts),
    dimnames = dimnames(train_counts)
  )

  adjusted_train[
    gene_keep_index,

  ] <- adjusted_fit_counts

  if (length(gene_remove_index) > 0L) {
    adjusted_train[
      gene_remove_index,

    ] <- train_counts[
      gene_remove_index,
      ,
      drop = FALSE
    ]
  }

  list(
    adjusted_train = adjusted_train,
    state = list(
      representation = "ComBatSeq_log1p",
      implementation = "training-estimated frozen ComBat-Seq",
      initial_genes = rownames(train_counts),
      fitted_genes = rownames(train_counts)[gene_keep_index],
      unchanged_genes = rownames(train_counts)[gene_remove_index],
      batch_levels = batch_levels,
      alpha_g = alpha_g,
      gamma_star = gamma_star,
      phi_hat = phi_hat,
      phi_star = phi_star
    )
  )
}

apply_frozen_combatseq <- function(new_counts,
                                   batch_new,
                                   state) {
  new_counts <- new_counts[
    state$initial_genes,
    ,
    drop = FALSE
  ]

  batch_new <- as.character(batch_new)

  unseen_batches <- setdiff(
    unique(batch_new),
    state$batch_levels
  )

  if (length(unseen_batches) > 0L) {
    stop(
      "Frozen ComBat-Seq cannot transform cohort(s) absent from fitting data: ",
      paste(unseen_batches, collapse = ", "),
      call. = FALSE
    )
  }

  fitted_gene_index <- match(
    state$fitted_genes,
    rownames(new_counts)
  )

  new_fit_counts <- new_counts[
    fitted_gene_index,
    ,
    drop = FALSE
  ]

  new_dge <- edgeR::DGEList(
    counts = new_fit_counts
  )

  new_library_offset <- edgeR::getOffset(
    new_dge
  )

  if (any(!is.finite(new_library_offset))) {
    stop(
      "Invalid frozen ComBat-Seq library-size offset in evaluation data.",
      call. = FALSE
    )
  }

  adjusted_fit_counts <- matrix(
    NA_real_,
    nrow = nrow(new_fit_counts),
    ncol = ncol(new_fit_counts),
    dimnames = dimnames(new_fit_counts)
  )

  for (batch_level in unique(batch_new)) {
    current_columns <- which(batch_new == batch_level)
    batch_index <- match(batch_level, state$batch_levels)

    base_mu <- exp(
      vec_to_matrix(
        state$alpha_g,
        length(current_columns)
      ) +
        matrix(
          rep(
            new_library_offset[current_columns],
            each = length(state$alpha_g)
          ),
          nrow = length(state$alpha_g),
          ncol = length(current_columns)
        )
    )

    old_mu <- base_mu *
      exp(
        vec_to_matrix(
          state$gamma_star[
            ,
            batch_index
          ],
          length(current_columns)
        )
      )

    adjusted_fit_counts[
      ,
      current_columns
    ] <- combatseq_match_quantiles(
      counts_sub = new_fit_counts[
        ,
        current_columns,
        drop = FALSE
      ],
      old_mu = old_mu,
      old_phi = state$phi_hat[
        ,
        batch_index
      ],
      new_mu = base_mu,
      new_phi = state$phi_star
    )
  }

  adjusted_new <- matrix(
    NA_real_,
    nrow = nrow(new_counts),
    ncol = ncol(new_counts),
    dimnames = dimnames(new_counts)
  )

  adjusted_new[
    fitted_gene_index,

  ] <- adjusted_fit_counts

  unchanged_gene_index <- setdiff(
    seq_len(nrow(new_counts)),
    fitted_gene_index
  )

  if (length(unchanged_gene_index) > 0L) {
    adjusted_new[
      unchanged_gene_index,

    ] <- new_counts[
      unchanged_gene_index,
      ,
      drop = FALSE
    ]
  }

  adjusted_new
}

fit_apply_frozen_combatseq_log1p <- function(train_counts,
                                             new_counts = NULL,
                                             train_metadata,
                                             new_metadata = NULL) {
  fit <- fit_frozen_combatseq(
    train_counts = train_counts,
    batch_train = train_metadata$cohort_group
  )

  train_log <- log1p(
    fit$adjusted_train
  )

  new_log <- NULL

  if (!is.null(new_counts)) {
    if (is.null(new_metadata)) {
      stop(
        "new_metadata is required for frozen ComBat-Seq projection.",
        call. = FALSE
      )
    }

    adjusted_new <- apply_frozen_combatseq(
      new_counts = new_counts,
      batch_new = new_metadata$cohort_group,
      state = fit$state
    )

    new_log <- log1p(
      adjusted_new
    )
  }

  if (anyNA(train_log) || !all(is.finite(train_log))) {
    stop(
      "Invalid frozen ComBat-Seq + log1p fitting matrix.",
      call. = FALSE
    )
  }

  if (
    !is.null(new_log) &&
      (anyNA(new_log) || !all(is.finite(new_log)))
  ) {
    stop(
      "Invalid frozen ComBat-Seq + log1p evaluation matrix.",
      call. = FALSE
    )
  }

  list(
    train = train_log,
    new = new_log,
    state = fit$state
  )
}


# ==============================================================================
# 10. SAMPLE-WISE LOG2(CPM + 1)
# ==============================================================================

fit_apply_log2_cpm_plus_1 <- function(train_counts,
                                      new_counts = NULL,
                                      min_total_count = MIN_TOTAL_TRAINING_COUNT,
                                      min_features = DEFAULT_TOP_N) {
  training_total_counts <- rowSums(train_counts)

  keep_gene <- (
    is.finite(training_total_counts) &
      training_total_counts >= min_total_count
  )

  if (sum(keep_gene) < min_features) {
    stop(
      "Fewer than ",
      min_features,
      " genes survived the fitting-only count filter for log2(CPM+1).",
      call. = FALSE
    )
  }

  genes <- rownames(train_counts)[keep_gene]

  transform_one <- function(counts) {
    counts <- counts[
      genes,
      ,
      drop = FALSE
    ]

    library_size <- colSums(counts)

    if (any(!is.finite(library_size)) || any(library_size <= 0)) {
      stop(
        "Invalid library size encountered during log2(CPM+1) transformation.",
        call. = FALSE
      )
    }

    cpm <- sweep(
      counts,
      MARGIN = 2L,
      STATS = library_size,
      FUN = "/"
    ) * 1e6

    log2(cpm + 1)
  }

  train <- transform_one(train_counts)
  new <- if (!is.null(new_counts)) transform_one(new_counts) else NULL

  list(
    train = train,
    new = new,
    state = list(
      representation = "log2_CPM_plus_1",
      genes = genes,
      min_total_count = min_total_count,
      transformation = "sample-wise log2(CPM + 1)"
    )
  )
}


# ==============================================================================
# 11. UNIFIED REPRESENTATION WRAPPER
# ==============================================================================

fit_apply_representation <- function(representation_id,
                                     train_counts,
                                     new_counts = NULL,
                                     train_metadata,
                                     new_metadata = NULL,
                                     seed = GLOBAL_SEED) {
  if (representation_id == "VST_only") {
    return(
      fit_apply_frozen_vst(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata
      )
    )
  }

  if (representation_id == "VST_fSVA") {
    return(
      fit_apply_vst_fsva(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        seed = seed
      )
    )
  }

  if (representation_id == "VST_ComBat") {
    return(
      fit_apply_vst_combat(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        new_metadata = new_metadata,
        reference_anchored = FALSE
      )
    )
  }

  if (representation_id == "VST_ComBatRef") {
    return(
      fit_apply_vst_combat(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        new_metadata = new_metadata,
        reference_anchored = TRUE
      )
    )
  }

  if (representation_id == "ComBatSeq_log1p") {
    return(
      fit_apply_frozen_combatseq_log1p(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        new_metadata = new_metadata
      )
    )
  }

  if (representation_id == "log2_CPM_plus_1") {
    return(
      fit_apply_log2_cpm_plus_1(
        train_counts = train_counts,
        new_counts = new_counts
      )
    )
  }

  stop(
    "Unknown representation identifier: ",
    representation_id,
    call. = FALSE
  )
}


# ==============================================================================
# 12. REPRESENTATION DIAGNOSTICS
# ==============================================================================

sanitize_for_pca <- function(mat) {
  mat <- coerce_numeric_matrix(mat)

  finite_gene <- rowSums(!is.finite(mat)) == 0L
  mat <- mat[finite_gene, , drop = FALSE]

  gene_var <- matrixStats::rowVars(mat)

  keep <- (
    is.finite(gene_var) &
      gene_var > sqrt(.Machine$double.eps)
  )

  mat[keep, , drop = FALSE]
}

run_pca_object <- function(mat) {
  mat <- sanitize_for_pca(mat)

  pr <- stats::prcomp(
    t(mat),
    center = TRUE,
    scale. = TRUE
  )

  list(
    scores = pr$x,
    variance = (pr$sdev^2) / sum(pr$sdev^2)
  )
}

silhouette_trajectory <- function(scores,
                                  labels,
                                  max_k = MAX_PCA_K) {
  labels <- droplevels(factor(labels))

  if (
    length(levels(labels)) < 2L ||
      nrow(scores) < 3L
  ) {
    return(data.table())
  }

  maximum_k <- min(
    max_k,
    ncol(scores),
    nrow(scores) - 1L
  )

  if (maximum_k < 2L) {
    return(data.table())
  }

  rbindlist(
    lapply(
      2:maximum_k,
      function(k) {
        xk <- scores[
          ,
          seq_len(k),
          drop = FALSE
        ]

        d <- stats::dist(xk)

        sil <- cluster::silhouette(
          as.integer(labels),
          d
        )

        data.table(
          k = k,
          mean_silhouette = mean(sil[, "sil_width"]),
          median_silhouette = median(sil[, "sil_width"])
        )
      }
    )
  )
}

choose_silhouette_k <- function(dt) {
  if (nrow(dt) == 0L) {
    return(
      data.table(
        k_star = NA_integer_,
        mean_silhouette = NA_real_,
        median_silhouette = NA_real_
      )
    )
  }

  max_mean <- max(dt$mean_silhouette, na.rm = TRUE)

  if (is.finite(max_mean) && max_mean > 0) {
    eligible <- dt[
      mean_silhouette >= SILHOUETTE_NEAR_MAX * max_mean
    ]

    chosen <- eligible[order(k)][1]
  } else {
    chosen <- dt[
      order(-mean_silhouette, k)
    ][1]
  }

  data.table(
    k_star = chosen$k,
    mean_silhouette = chosen$mean_silhouette,
    median_silhouette = chosen$median_silhouette
  )
}

calculate_depth_dependence <- function(representations,
                                       raw_counts,
                                       metadata) {
  overall_rows <- list()
  cohort_rows <- list()

  for (representation_id in names(representations)) {
    mat <- representations[[representation_id]]

    matched_raw <- raw_counts[
      rownames(mat),
      colnames(mat),
      drop = FALSE
    ]

    raw_library_size <- colSums(matched_raw)
    transformed_total_signal <- colSums(mat)

    sample_dt <- data.table(
      sample_id = colnames(mat),
      raw_library_size = raw_library_size,
      transformed_total_signal = transformed_total_signal
    )

    meta_dt <- data.table(
      sample_id = metadata$sample_id,
      cohort = as.character(metadata$cohort_group)
    )

    sample_dt <- merge(
      sample_dt,
      meta_dt,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )

    valid <- (
      is.finite(sample_dt$raw_library_size) &
        is.finite(sample_dt$transformed_total_signal)
    )

    pooled_test <- if (
      sum(valid) >= 3L &&
        length(unique(sample_dt$raw_library_size[valid])) > 1L &&
        length(unique(sample_dt$transformed_total_signal[valid])) > 1L
    ) {
      suppressWarnings(
        stats::cor.test(
          sample_dt$raw_library_size[valid],
          sample_dt$transformed_total_signal[valid],
          method = "spearman",
          exact = FALSE
        )
      )
    } else {
      NULL
    }

    overall_rows[[representation_id]] <- data.table(
      representation = representation_id,
      n = sum(valid),
      rho = if (!is.null(pooled_test)) unname(pooled_test$estimate) else NA_real_,
      p_value = if (!is.null(pooled_test)) pooled_test$p.value else NA_real_
    )

    cohort_rows_current <- rbindlist(
      lapply(
        unique(sample_dt$cohort),
        function(cohort_name) {
          current <- sample_dt[cohort == cohort_name]

          valid_current <- (
            is.finite(current$raw_library_size) &
              is.finite(current$transformed_total_signal)
          )

          current_test <- if (
            sum(valid_current) >= 3L &&
              length(unique(current$raw_library_size[valid_current])) > 1L &&
              length(unique(current$transformed_total_signal[valid_current])) > 1L
          ) {
            suppressWarnings(
              stats::cor.test(
                current$raw_library_size[valid_current],
                current$transformed_total_signal[valid_current],
                method = "spearman",
                exact = FALSE
              )
            )
          } else {
            NULL
          }

          data.table(
            representation = representation_id,
            cohort = cohort_name,
            n = sum(valid_current),
            rho = if (!is.null(current_test)) unname(current_test$estimate) else NA_real_,
            p_value = if (!is.null(current_test)) current_test$p.value else NA_real_
          )
        }
      ),
      fill = TRUE
    )

    cohort_rows[[representation_id]] <- cohort_rows_current
  }

  overall <- rbindlist(overall_rows, fill = TRUE)
  by_cohort <- rbindlist(cohort_rows, fill = TRUE)

  overall[, p_adjusted := stats::p.adjust(p_value, method = "BH")]
  by_cohort[, p_adjusted := stats::p.adjust(p_value, method = "BH")]

  list(
    overall = overall,
    by_cohort = by_cohort
  )
}

run_representation_diagnostics <- function(representations,
                                           raw_counts,
                                           metadata,
                                           output_dir,
                                           prefix) {
  make_dir(output_dir)

  pca_rows <- list()
  variance_rows <- list()
  silhouette_rows <- list()
  silhouette_summary_rows <- list()

  for (representation_index in seq_along(representations)) {
    representation_id <- names(representations)[[representation_index]]
    mat <- representations[[representation_id]]

    pca <- run_pca_object(mat)

    scores <- pca$scores[
      metadata$sample_id,
      ,
      drop = FALSE
    ]

    pca_rows[[representation_id]] <- data.table(
      sample_id = rownames(scores),
      PC1 = scores[, 1],
      PC2 = scores[, 2],
      representation = representation_id,
      cohort = as.character(metadata$cohort_group),
      outcome = as.character(metadata$outcome)
    )

    variance_rows[[representation_id]] <- data.table(
      representation = representation_id,
      PC = seq_along(pca$variance),
      variance_fraction = pca$variance
    )

    cohort_sil <- silhouette_trajectory(
      scores,
      metadata$cohort_group
    )

    outcome_sil <- silhouette_trajectory(
      scores,
      metadata$outcome
    )

    if (nrow(cohort_sil) > 0L) {
      cohort_sil[, `:=`(
        representation = representation_id,
        grouping = "cohort"
      )]
    }

    if (nrow(outcome_sil) > 0L) {
      outcome_sil[, `:=`(
        representation = representation_id,
        grouping = "outcome"
      )]
    }

    current_sil <- rbind(
      cohort_sil,
      outcome_sil,
      fill = TRUE
    )

    silhouette_rows[[representation_id]] <- current_sil

    if (nrow(cohort_sil) > 0L) {
      cohort_star <- choose_silhouette_k(cohort_sil)
      cohort_star[, `:=`(
        representation = representation_id,
        grouping = "cohort"
      )]

      silhouette_summary_rows[[paste0(representation_id, "_cohort")]] <- cohort_star
    }

    if (nrow(outcome_sil) > 0L) {
      outcome_star <- choose_silhouette_k(outcome_sil)
      outcome_star[, `:=`(
        representation = representation_id,
        grouping = "outcome"
      )]

      silhouette_summary_rows[[paste0(representation_id, "_outcome")]] <- outcome_star
    }
  }

  pca_dt <- rbindlist(pca_rows, fill = TRUE)
  variance_dt <- rbindlist(variance_rows, fill = TRUE)
  silhouette_dt <- rbindlist(silhouette_rows, fill = TRUE)
  silhouette_summary <- rbindlist(silhouette_summary_rows, fill = TRUE)

  depth <- calculate_depth_dependence(
    representations = representations,
    raw_counts = raw_counts,
    metadata = metadata
  )

  safe_write_csv(
    pca_dt,
    file.path(output_dir, paste0(prefix, "_pca_scores.csv"))
  )

  safe_write_csv(
    variance_dt,
    file.path(output_dir, paste0(prefix, "_pca_variance.csv"))
  )

  safe_write_csv(
    silhouette_dt,
    file.path(output_dir, paste0(prefix, "_silhouette_all_k.csv"))
  )

  safe_write_csv(
    silhouette_summary,
    file.path(output_dir, paste0(prefix, "_silhouette_summary.csv"))
  )

  safe_write_csv(
    depth$overall,
    file.path(output_dir, paste0(prefix, "_depth_dependence_overall.csv"))
  )

  safe_write_csv(
    depth$by_cohort,
    file.path(output_dir, paste0(prefix, "_depth_dependence_by_cohort.csv"))
  )

  p_pca_cohort <- ggplot(
    pca_dt,
    aes(
      x = PC1,
      y = PC2,
      colour = cohort
    )
  ) +
    geom_point(alpha = 0.75, size = 1.4) +
    facet_wrap(~ representation, scales = "free") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "PC1",
      y = "PC2",
      colour = "Cohort"
    )

  p_pca_outcome <- ggplot(
    pca_dt,
    aes(
      x = PC1,
      y = PC2,
      colour = outcome
    )
  ) +
    geom_point(alpha = 0.75, size = 1.4) +
    facet_wrap(~ representation, scales = "free") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "PC1",
      y = "PC2",
      colour = "Outcome"
    )

  ggsave(
    file.path(output_dir, paste0(prefix, "_pca_by_cohort.png")),
    p_pca_cohort,
    width = 180,
    height = 130,
    units = "mm",
    dpi = 600
  )

  ggsave(
    file.path(output_dir, paste0(prefix, "_pca_by_outcome.png")),
    p_pca_outcome,
    width = 180,
    height = 130,
    units = "mm",
    dpi = 600
  )

  if (nrow(silhouette_dt) > 0L) {
    p_sil <- ggplot(
      silhouette_dt,
      aes(
        x = k,
        y = mean_silhouette,
        colour = grouping
      )
    ) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ representation) +
      theme_bw(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold")
      ) +
      labs(
        x = "Number of principal components",
        y = "Mean silhouette width",
        colour = "Grouping"
      )

    ggsave(
      file.path(output_dir, paste0(prefix, "_silhouette_profiles.png")),
      p_sil,
      width = 180,
      height = 125,
      units = "mm",
      dpi = 600
    )
  }

  invisible(
    list(
      pca = pca_dt,
      variance = variance_dt,
      silhouette = silhouette_dt,
      silhouette_summary = silhouette_summary,
      depth = depth
    )
  )
}


# ==============================================================================
# 13. SAVE A REPRESENTATION BUNDLE
# ==============================================================================

save_representation_bundle <- function(fit,
                                       output_dir,
                                       prefix,
                                       train_metadata,
                                       new_metadata = NULL) {
  make_dir(output_dir)

  safe_save_rds(
    fit$train,
    file.path(output_dir, paste0(prefix, "_development_matrix.rds"))
  )

  if (!is.null(fit$new)) {
    safe_save_rds(
      fit$new,
      file.path(output_dir, paste0(prefix, "_evaluation_matrix.rds"))
    )
  }

  safe_save_rds(
    fit$state,
    file.path(output_dir, paste0(prefix, "_fitted_state.rds"))
  )

  safe_write_csv(
    train_metadata,
    file.path(output_dir, paste0(prefix, "_development_metadata.csv"))
  )

  if (!is.null(new_metadata)) {
    safe_write_csv(
      new_metadata,
      file.path(output_dir, paste0(prefix, "_evaluation_metadata.csv"))
    )
  }

  summary <- data.table(
    representation = if (!is.null(fit$state$representation)) fit$state$representation else prefix,
    development_n = ncol(fit$train),
    development_genes = nrow(fit$train),
    evaluation_n = if (!is.null(fit$new)) ncol(fit$new) else 0L,
    evaluation_genes = if (!is.null(fit$new)) nrow(fit$new) else NA_integer_
  )

  safe_write_csv(
    summary,
    file.path(output_dir, paste0(prefix, "_bundle_summary.csv"))
  )

  invisible(summary)
}


# ==============================================================================
# 14. FULL-DEVELOPMENT REPRESENTATION SPECIFICATION
# ==============================================================================

run_full_development_representation_stage <- function(disease_object,
                                                      split_object) {
  disease_id <- disease_object$disease_id

  output_dir <- file.path(
    STAGE_ROOT,
    "01_full_development",
    disease_id
  )

  make_dir(output_dir)

  message_section(
    paste0(
      "Full-development representation specification | ",
      disease_id
    )
  )

  development_counts <- split_object$development_counts
  development_metadata <- split_object$development_metadata

  vst_only <- fit_apply_frozen_vst(
    train_counts = development_counts,
    new_counts = NULL,
    train_metadata = development_metadata
  )

  vst_combat <- fit_apply_vst_combat(
    train_counts = development_counts,
    new_counts = NULL,
    train_metadata = development_metadata,
    new_metadata = NULL,
    reference_anchored = FALSE
  )

  vst_combat_ref <- fit_apply_vst_combat(
    train_counts = development_counts,
    new_counts = NULL,
    train_metadata = development_metadata,
    new_metadata = NULL,
    reference_anchored = TRUE
  )

  vst_fsva <- fit_apply_vst_fsva(
    train_counts = development_counts,
    new_counts = NULL,
    train_metadata = development_metadata,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      disease_id,
      "full_development_fsva"
    )
  )

  combatseq <- fit_standard_combatseq_diagnostic(
    train_counts = development_counts,
    train_metadata = development_metadata
  )

  raw_genes <- vst_only$state$genes

  representations <- list(
    RAW_counts = development_counts[
      raw_genes,
      development_metadata$sample_id,
      drop = FALSE
    ],
    VST_only = vst_only$train,
    VST_ComBat = vst_combat$train,
    VST_ComBatRef = vst_combat_ref$train,
    VST_fSVA = vst_fsva$train,
    ComBatSeq_log1p = combatseq$train
  )

  safe_save_rds(
    representations,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_development_representation_matrices.rds"
      )
    )
  )

  states <- list(
    VST_only = vst_only$state,
    VST_ComBat = vst_combat$state,
    VST_ComBatRef = vst_combat_ref$state,
    VST_fSVA = vst_fsva$state,
    ComBatSeq_log1p = combatseq$state
  )

  safe_save_rds(
    states,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_development_representation_states.rds"
      )
    )
  )

  diagnostic_raw_counts <- development_counts[
    ,
    development_metadata$sample_id,
    drop = FALSE
  ]

  run_representation_diagnostics(
    representations = representations,
    raw_counts = diagnostic_raw_counts,
    metadata = development_metadata,
    output_dir = file.path(output_dir, "diagnostics"),
    prefix = disease_id
  )

  state_summary <- data.table(
    representation = names(representations),
    genes = vapply(representations, nrow, integer(1)),
    samples = vapply(representations, ncol, integer(1))
  )

  safe_write_csv(
    state_summary,
    file.path(
      output_dir,
      paste0(
        disease_id,
        "_representation_summary.csv"
      )
    )
  )

  invisible(
    list(
      representations = representations,
      states = states
    )
  )
}


# ==============================================================================
# 15. EXPERIMENT 1: SAMPLE-SUPPORT PERTURBATION
# ==============================================================================

run_experiment_1 <- function(disease_objects) {
  experiment_dir <- file.path(
    STAGE_ROOT,
    "02_experiment_1_sample_support"
  )

  make_dir(experiment_dir)

  experiment_diseases <- intersect(
    c("TB", "RA", "UC", "CD"),
    names(disease_objects)
  )

  if (length(experiment_diseases) == 0L) {
    warning(
      "Experiment 1 skipped because none of TB/RA/UC/CD were loaded.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  exclusion_table <- data.table(
    disease = "T1D",
    included = FALSE,
    reason = "Available T1D sample support is below the common 82-observation design."
  )

  safe_write_csv(
    exclusion_table,
    file.path(experiment_dir, "T1D_exclusion.csv")
  )

  for (disease_id in experiment_diseases) {
    disease_object <- disease_objects[[disease_id]]
    counts <- disease_object$counts
    metadata <- disease_object$metadata

    output_dir <- file.path(experiment_dir, disease_id)
    split_dir <- file.path(output_dir, "00_manifests")
    matrix_dir <- file.path(output_dir, "01_representations")

    make_dir(split_dir)
    make_dir(matrix_dir)

    message_section(
      paste0(
        "Experiment 1 | equal-size sample support | ",
        disease_id
      )
    )

    selected_cases <- select_ids_by_outcome_and_cohort(
      metadata = metadata,
      outcome_value = POSITIVE_CLASS,
      target = 41L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "equal_size",
        "case"
      ),
      retain_at_least_one = FALSE
    )

    selected_controls <- select_ids_by_outcome_and_cohort(
      metadata = metadata,
      outcome_value = NEGATIVE_CLASS,
      target = 41L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "equal_size",
        "control"
      ),
      retain_at_least_one = FALSE
    )

    equal_ids <- c(
      selected_cases$ids,
      selected_controls$ids
    )

    equal_metadata <- metadata[
      match(equal_ids, metadata$sample_id),
      ,
      drop = FALSE
    ]

    equal_counts <- counts[
      ,
      equal_metadata$sample_id,
      drop = FALSE
    ]

    if (
      sum(equal_metadata$outcome == POSITIVE_CLASS) != 41L ||
        sum(equal_metadata$outcome == NEGATIVE_CLASS) != 41L ||
        nrow(equal_metadata) != 82L ||
        anyDuplicated(equal_metadata$sample_id) > 0L
    ) {
      stop(
        disease_id,
        " equal-size 41/41 selection audit failed.",
        call. = FALSE
      )
    }

    equal_allocation <- rbind(
      selected_cases$allocation,
      selected_controls$allocation,
      fill = TRUE
    )

    safe_write_csv(
      equal_allocation,
      file.path(
        split_dir,
        paste0(
          disease_id,
          "_equal_size_41v41_allocation.csv"
        )
      )
    )

    heldout_cases <- select_ids_by_outcome_and_cohort(
      metadata = equal_metadata,
      outcome_value = POSITIVE_CLASS,
      target = 8L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "equal_size",
        "held_out_case"
      ),
      retain_at_least_one = TRUE
    )

    heldout_controls <- select_ids_by_outcome_and_cohort(
      metadata = equal_metadata,
      outcome_value = NEGATIVE_CLASS,
      target = 8L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "equal_size",
        "held_out_control"
      ),
      retain_at_least_one = TRUE
    )

    heldout_ids <- c(
      heldout_cases$ids,
      heldout_controls$ids
    )

    development_ids <- setdiff(
      equal_metadata$sample_id,
      heldout_ids
    )

    development_metadata <- equal_metadata[
      match(development_ids, equal_metadata$sample_id),
      ,
      drop = FALSE
    ]

    heldout_metadata <- equal_metadata[
      match(heldout_ids, equal_metadata$sample_id),
      ,
      drop = FALSE
    ]

    development_counts <- equal_counts[
      ,
      development_metadata$sample_id,
      drop = FALSE
    ]

    heldout_counts <- equal_counts[
      ,
      heldout_metadata$sample_id,
      drop = FALSE
    ]

    assert_disjoint_samples(
      development_metadata$sample_id,
      heldout_metadata$sample_id
    )

    if (
      nrow(development_metadata) != 66L ||
        sum(development_metadata$outcome == POSITIVE_CLASS) != 33L ||
        sum(development_metadata$outcome == NEGATIVE_CLASS) != 33L ||
        nrow(heldout_metadata) != 16L ||
        sum(heldout_metadata$outcome == POSITIVE_CLASS) != 8L ||
        sum(heldout_metadata$outcome == NEGATIVE_CLASS) != 8L
    ) {
      stop(
        disease_id,
        " equal-size 66/16 partition audit failed.",
        call. = FALSE
      )
    }

    manifest <- rbind(
      data.table(
        disease = disease_id,
        sample_id = development_metadata$sample_id,
        cohort = as.character(development_metadata$cohort_group),
        outcome = as.character(development_metadata$outcome),
        role = "development"
      ),
      data.table(
        disease = disease_id,
        sample_id = heldout_metadata$sample_id,
        cohort = as.character(heldout_metadata$cohort_group),
        outcome = as.character(heldout_metadata$outcome),
        role = "held_out"
      )
    )

    safe_write_csv(
      manifest,
      file.path(
        split_dir,
        paste0(
          disease_id,
          "_equal_size_82_manifest.csv"
        )
      )
    )

    representation_ids <- c(
      "VST_fSVA",
      "VST_ComBatRef",
      "ComBatSeq_log1p"
    )

    development_representations <- list()

    for (representation_index in seq_along(representation_ids)) {
      representation_id <- representation_ids[[representation_index]]

      fit <- fit_apply_representation(
        representation_id = representation_id,
        train_counts = development_counts,
        new_counts = heldout_counts,
        train_metadata = development_metadata,
        new_metadata = heldout_metadata,
        seed = make_trace_subseed(
          GLOBAL_SEED,
          disease_id,
          "experiment1",
          representation_id
        )
      )

      representation_output_dir <- file.path(
        matrix_dir,
        representation_id
      )

      save_representation_bundle(
        fit = fit,
        output_dir = representation_output_dir,
        prefix = paste0(disease_id, "_", representation_id),
        train_metadata = development_metadata,
        new_metadata = heldout_metadata
      )

      development_representations[[representation_id]] <- fit$train
    }

    diagnostics_raw <- development_counts[
      ,
      development_metadata$sample_id,
      drop = FALSE
    ]

    run_representation_diagnostics(
      representations = development_representations,
      raw_counts = diagnostics_raw,
      metadata = development_metadata,
      output_dir = file.path(output_dir, "02_diagnostics"),
      prefix = paste0(disease_id, "_equal_size")
    )

    feasibility <- data.table(
      representation = representation_ids,
      ordinary_cv_with_known_cohort_levels = c(TRUE, TRUE, TRUE),
      frozen_heldout_with_known_cohort_levels = c(TRUE, TRUE, TRUE),
      strict_unseen_cohort_loco = c(TRUE, FALSE, FALSE),
      reason = c(
        "fSVA projects a development-estimated latent structure to new observations.",
        "ComBat-Ref projection requires target cohort levels to occur in fitting data.",
        "Frozen ComBat-Seq projection requires target cohort levels to occur in fitting data."
      )
    )

    safe_write_csv(
      feasibility,
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_representation_feasibility.csv"
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 16. CROSS-DISEASE POOLING HELPERS
# ==============================================================================

prepare_cross_disease_component <- function(disease_id,
                                            split_object,
                                            common_genes,
                                            role = "development") {
  if (role == "development") {
    counts <- split_object$development_counts
    metadata <- split_object$development_metadata
  } else if (role == "held_out") {
    counts <- split_object$held_out_counts
    metadata <- split_object$held_out_metadata
  } else {
    stop("Unknown cross-disease role.", call. = FALSE)
  }

  counts <- counts[
    common_genes,
    metadata$sample_id,
    drop = FALSE
  ]

  metadata <- metadata
  metadata$original_sample_id <- metadata$sample_id
  metadata$sample_id <- paste(
    disease_id,
    metadata$original_sample_id,
    sep = "__"
  )

  colnames(counts) <- metadata$sample_id
  metadata$disease <- disease_id
  metadata$cohort_group <- factor(
    paste(
      disease_id,
      as.character(metadata$cohort_group),
      sep = "__"
    )
  )

  list(
    counts = counts,
    metadata = metadata
  )
}

pool_cross_disease_components <- function(components) {
  genes <- rownames(components[[1]]$counts)

  for (component in components) {
    if (!identical(rownames(component$counts), genes)) {
      stop(
        "Cross-disease components have inconsistent feature order.",
        call. = FALSE
      )
    }
  }

  counts <- do.call(
    cbind,
    lapply(
      components,
      function(x) x$counts
    )
  )

  metadata <- do.call(
    rbind,
    lapply(
      components,
      function(x) x$metadata
    )
  )

  metadata <- as.data.frame(
    metadata,
    check.names = FALSE
  )

  metadata <- metadata[
    match(colnames(counts), metadata$sample_id),
    ,
    drop = FALSE
  ]

  if (!identical(colnames(counts), metadata$sample_id)) {
    stop(
      "Cross-disease pooled counts and metadata are misaligned.",
      call. = FALSE
    )
  }

  list(
    counts = counts,
    metadata = metadata
  )
}


# ==============================================================================
# 17. EXPERIMENT 2: DISEASE-COMPOSITION / LODO
# ==============================================================================

run_experiment_2 <- function(disease_objects,
                             general_splits) {
  required <- c("TB", "RA", "UC", "CD", "T1D")

  if (!all(required %in% names(disease_objects))) {
    stop(
      "Experiment 2 requires TB, RA, UC, CD and T1D inputs.",
      call. = FALSE
    )
  }

  output_dir <- file.path(
    STAGE_ROOT,
    "03_experiment_2_disease_composition"
  )

  make_dir(output_dir)

  message_section("Experiment 2 | disease composition and LODO")

  development_gene_sets <- lapply(
    required,
    function(disease_id) {
      rownames(
        general_splits[[disease_id]]$development_counts
      )
    }
  )

  names(development_gene_sets) <- required

  common_genes <- Reduce(
    intersect,
    development_gene_sets
  )

  common_genes <- sort(common_genes)

  safe_write_csv(
    data.table(
      rank = seq_along(common_genes),
      gene = common_genes
    ),
    file.path(
      output_dir,
      "cross_disease_common_genes.csv"
    )
  )

  gene_audit <- data.table(
    disease = c(required, "ALL5"),
    available_genes = c(
      vapply(
        development_gene_sets,
        length,
        integer(1)
      ),
      length(common_genes)
    )
  )

  safe_write_csv(
    gene_audit,
    file.path(
      output_dir,
      "cross_disease_gene_space_audit.csv"
    )
  )

  audit_expected_value(
    observed = length(common_genes),
    expected = EXPECTED_COMMON_CROSS_DISEASE_GENES,
    label = "Cross-disease common feature space"
  )

  development_components <- lapply(
    required,
    function(disease_id) {
      prepare_cross_disease_component(
        disease_id = disease_id,
        split_object = general_splits[[disease_id]],
        common_genes = common_genes,
        role = "development"
      )
    }
  )

  names(development_components) <- required

  tb_heldout_component <- prepare_cross_disease_component(
    disease_id = "TB",
    split_object = general_splits[["TB"]],
    common_genes = common_genes,
    role = "held_out"
  )

  model_a_diseases <- c("RA", "UC", "CD", "T1D")
  model_b_diseases <- required

  model_a <- pool_cross_disease_components(
    development_components[model_a_diseases]
  )

  model_b <- pool_cross_disease_components(
    development_components[model_b_diseases]
  )

  audit_expected_value(
    observed = ncol(model_a$counts),
    expected = EXPECTED_MODEL_A_N,
    label = "Cross-disease Model A development n"
  )

  audit_expected_value(
    observed = ncol(model_b$counts),
    expected = EXPECTED_MODEL_B_N,
    label = "Cross-disease Model B development n"
  )

  model_a_fit <- fit_apply_vst_fsva(
    train_counts = model_a$counts,
    new_counts = tb_heldout_component$counts,
    train_metadata = model_a$metadata,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "experiment2",
      "model_a"
    ),
    include_cohort_in_vst_design = FALSE
  )

  model_b_fit <- fit_apply_vst_fsva(
    train_counts = model_b$counts,
    new_counts = tb_heldout_component$counts,
    train_metadata = model_b$metadata,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      "experiment2",
      "model_b"
    ),
    include_cohort_in_vst_design = FALSE
  )

  save_representation_bundle(
    fit = model_a_fit,
    output_dir = file.path(output_dir, "Model_A_TB_naive"),
    prefix = "Model_A_VST_fSVA",
    train_metadata = model_a$metadata,
    new_metadata = tb_heldout_component$metadata
  )

  save_representation_bundle(
    fit = model_b_fit,
    output_dir = file.path(output_dir, "Model_B_ALL5"),
    prefix = "Model_B_VST_fSVA",
    train_metadata = model_b$metadata,
    new_metadata = tb_heldout_component$metadata
  )

  model_composition <- data.table(
    model = c("A", "B"),
    diseases = c(
      paste(model_a_diseases, collapse = " + "),
      paste(model_b_diseases, collapse = " + ")
    ),
    development_n = c(
      ncol(model_a$counts),
      ncol(model_b$counts)
    ),
    common_genes = length(common_genes),
    target = "Frozen TB held-out population"
  )

  safe_write_csv(
    model_composition,
    file.path(
      output_dir,
      "model_a_model_b_composition.csv"
    )
  )

  lodo_summary_rows <- list()

  for (target_index in seq_along(required)) {
    target_disease <- required[[target_index]]
    source_diseases <- setdiff(required, target_disease)

    source <- pool_cross_disease_components(
      development_components[source_diseases]
    )

    target <- development_components[[target_disease]]

    fit <- fit_apply_vst_fsva(
      train_counts = source$counts,
      new_counts = target$counts,
      train_metadata = source$metadata,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        "experiment2",
        "LODO",
        target_disease
      ),
      include_cohort_in_vst_design = FALSE
    )

    target_dir <- file.path(
      output_dir,
      "LODO",
      paste0("target_", target_disease)
    )

    save_representation_bundle(
      fit = fit,
      output_dir = target_dir,
      prefix = paste0(
        "LODO_",
        target_disease,
        "_VST_fSVA"
      ),
      train_metadata = source$metadata,
      new_metadata = target$metadata
    )

    lodo_summary_rows[[target_disease]] <- data.table(
      target_disease = target_disease,
      source_diseases = paste(source_diseases, collapse = " + "),
      source_n = ncol(source$counts),
      target_n = ncol(target$counts),
      genes = nrow(fit$train),
      representation = "VST + fSVA"
    )
  }

  safe_write_csv(
    rbindlist(lodo_summary_rows, fill = TRUE),
    file.path(
      output_dir,
      "LODO_representation_summary.csv"
    )
  )

  invisible(TRUE)
}


# ==============================================================================
# 18. EXPERIMENT 3: COHORT-COMPOSITION / LOCO
# ==============================================================================

run_experiment_3 <- function(disease_objects,
                             general_splits) {
  experiment_dir <- file.path(
    STAGE_ROOT,
    "04_experiment_3_cohort_composition"
  )

  make_dir(experiment_dir)

  for (disease_id in names(disease_objects)) {
    split <- general_splits[[disease_id]]

    output_dir <- file.path(
      experiment_dir,
      disease_id
    )

    make_dir(output_dir)

    message_section(
      paste0(
        "Experiment 3 | cohort composition | ",
        disease_id
      )
    )

    pooled_vst <- fit_apply_frozen_vst(
      train_counts = split$development_counts,
      new_counts = split$held_out_counts,
      train_metadata = split$development_metadata
    )

    pooled_logcpm <- fit_apply_log2_cpm_plus_1(
      train_counts = split$development_counts,
      new_counts = split$held_out_counts
    )

    save_representation_bundle(
      fit = pooled_vst,
      output_dir = file.path(output_dir, "frozen_heldout", "VST_only"),
      prefix = paste0(disease_id, "_VST_only"),
      train_metadata = split$development_metadata,
      new_metadata = split$held_out_metadata
    )

    save_representation_bundle(
      fit = pooled_logcpm,
      output_dir = file.path(output_dir, "frozen_heldout", "log2_CPM_plus_1"),
      prefix = paste0(disease_id, "_log2_CPM_plus_1"),
      train_metadata = split$development_metadata,
      new_metadata = split$held_out_metadata
    )

    target_cohorts <- unique(
      as.character(
        split$development_metadata$cohort_group
      )
    )

    loco_rows <- list()

    for (cohort_index in seq_along(target_cohorts)) {
      target_cohort <- target_cohorts[[cohort_index]]

      target_rows <- (
        as.character(
          split$development_metadata$cohort_group
        ) == target_cohort
      )

      source_metadata <- split$development_metadata[
        !target_rows,
        ,
        drop = FALSE
      ]

      target_metadata <- split$development_metadata[
        target_rows,
        ,
        drop = FALSE
      ]

      source_counts <- split$development_counts[
        ,
        source_metadata$sample_id,
        drop = FALSE
      ]

      target_counts <- split$development_counts[
        ,
        target_metadata$sample_id,
        drop = FALSE
      ]

      assert_disjoint_samples(
        source_metadata$sample_id,
        target_metadata$sample_id,
        development_label = "LOCO source",
        evaluation_label = "LOCO target"
      )

      vst_fit <- fit_apply_frozen_vst(
        train_counts = source_counts,
        new_counts = target_counts,
        train_metadata = source_metadata
      )

      logcpm_fit <- fit_apply_log2_cpm_plus_1(
        train_counts = source_counts,
        new_counts = target_counts
      )

      cohort_safe <- safe_file_name(target_cohort)

      save_representation_bundle(
        fit = vst_fit,
        output_dir = file.path(
          output_dir,
          "LOCO",
          cohort_safe,
          "VST_only"
        ),
        prefix = paste0(
          disease_id,
          "_LOCO_",
          cohort_safe,
          "_VST_only"
        ),
        train_metadata = source_metadata,
        new_metadata = target_metadata
      )

      save_representation_bundle(
        fit = logcpm_fit,
        output_dir = file.path(
          output_dir,
          "LOCO",
          cohort_safe,
          "log2_CPM_plus_1"
        ),
        prefix = paste0(
          disease_id,
          "_LOCO_",
          cohort_safe,
          "_log2_CPM_plus_1"
        ),
        train_metadata = source_metadata,
        new_metadata = target_metadata
      )

      loco_rows[[target_cohort]] <- data.table(
        disease = disease_id,
        target_cohort = target_cohort,
        source_n = nrow(source_metadata),
        target_n = nrow(target_metadata),
        source_cohorts = length(unique(source_metadata$cohort_group)),
        vst_genes = nrow(vst_fit$train),
        log2cpm_genes = nrow(logcpm_fit$train)
      )
    }

    safe_write_csv(
      rbindlist(loco_rows, fill = TRUE),
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_LOCO_representation_summary.csv"
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 19. EXPERIMENT 4: MODEL-ARCHITECTURE PERTURBATION
# ==============================================================================

run_experiment_4 <- function(disease_objects,
                             general_splits) {
  experiment_dir <- file.path(
    STAGE_ROOT,
    "05_experiment_4_model_architecture"
  )

  make_dir(experiment_dir)

  for (disease_id in names(disease_objects)) {
    split <- general_splits[[disease_id]]

    output_dir <- file.path(
      experiment_dir,
      disease_id
    )

    make_dir(output_dir)

    message_section(
      paste0(
        "Experiment 4 | model architecture | ",
        disease_id
      )
    )

    matched_representation <- fit_apply_log2_cpm_plus_1(
      train_counts = split$development_counts,
      new_counts = split$held_out_counts
    )

    save_representation_bundle(
      fit = matched_representation,
      output_dir = file.path(output_dir, "frozen_diagnostic_test"),
      prefix = paste0(
        disease_id,
        "_architecture_log2_CPM_plus_1"
      ),
      train_metadata = split$development_metadata,
      new_metadata = split$held_out_metadata
    )

    outer_fold_map <- make_stratified_fold_map(
      metadata = split$development_metadata,
      folds = 5L,
      seed = make_trace_subseed(
        GLOBAL_SEED,
        disease_id,
        "architecture_outer_cv"
      )
    )

    safe_write_csv(
      outer_fold_map,
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_architecture_outer_fold_map.csv"
        )
      )
    )

    matched_input_audit <- data.table(
      disease = disease_id,
      representation = "log2(CPM + 1)",
      development_n = ncol(matched_representation$train),
      frozen_test_n = ncol(matched_representation$new),
      available_genes_before_feature_selection = nrow(matched_representation$train),
      downstream_selected_genes = 50L,
      architecture_1 = "LightGBM",
      architecture_2 = "TabPFN",
      identical_sample_population = TRUE,
      identical_feature_order_required = TRUE
    )

    safe_write_csv(
      matched_input_audit,
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_matched_architecture_input_audit.csv"
        )
      )
    )

    target_cohorts <- unique(
      as.character(
        split$development_metadata$cohort_group
      )
    )

    loco_rows <- list()

    for (cohort_index in seq_along(target_cohorts)) {
      target_cohort <- target_cohorts[[cohort_index]]

      target_rows <- (
        as.character(
          split$development_metadata$cohort_group
        ) == target_cohort
      )

      source_metadata <- split$development_metadata[
        !target_rows,
        ,
        drop = FALSE
      ]

      target_metadata <- split$development_metadata[
        target_rows,
        ,
        drop = FALSE
      ]

      source_counts <- split$development_counts[
        ,
        source_metadata$sample_id,
        drop = FALSE
      ]

      target_counts <- split$development_counts[
        ,
        target_metadata$sample_id,
        drop = FALSE
      ]

      fit <- fit_apply_log2_cpm_plus_1(
        train_counts = source_counts,
        new_counts = target_counts
      )

      cohort_safe <- safe_file_name(target_cohort)

      save_representation_bundle(
        fit = fit,
        output_dir = file.path(
          output_dir,
          "LOCO",
          cohort_safe
        ),
        prefix = paste0(
          disease_id,
          "_architecture_LOCO_",
          cohort_safe
        ),
        train_metadata = source_metadata,
        new_metadata = target_metadata
      )

      loco_rows[[target_cohort]] <- data.table(
        disease = disease_id,
        target_cohort = target_cohort,
        source_n = nrow(source_metadata),
        target_n = nrow(target_metadata),
        genes_before_RF_feature_selection = nrow(fit$train),
        downstream_selected_genes = 50L,
        representation = "log2(CPM + 1)"
      )
    }

    safe_write_csv(
      rbindlist(loco_rows, fill = TRUE),
      file.path(
        output_dir,
        paste0(
          disease_id,
          "_architecture_LOCO_representation_summary.csv"
        )
      )
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 20. LOAD ALL REQUESTED DISEASES AND CREATE FROZEN GENERAL MANIFESTS
# ==============================================================================

disease_objects <- list()
general_manifests <- list()
general_splits <- list()

for (disease_id in REQUESTED_DISEASES) {
  disease_objects[[disease_id]] <- load_disease_data(disease_id)

  manifest_path <- file.path(
    disease_objects[[disease_id]]$manifest_dir,
    paste0(
      disease_id,
      "_development_heldout_manifest.csv"
    )
  )

  development_fraction <- as.numeric(
    get_config_value(
      config,
      "partitioning.development_fraction",
      0.80
    )
  )

  manifest <- create_general_development_manifest(
    metadata = disease_objects[[disease_id]]$metadata,
    disease_id = disease_id,
    manifest_path = manifest_path,
    development_fraction = development_fraction,
    seed = make_trace_subseed(
      GLOBAL_SEED,
      disease_id,
      "general_development_split"
    ),
    overwrite = OVERWRITE_MANIFESTS
  )

  general_manifests[[disease_id]] <- manifest

  general_splits[[disease_id]] <- subset_from_manifest(
    disease_objects[[disease_id]],
    manifest
  )

  split_audit <- data.table(
    disease = disease_id,
    development_n = nrow(general_splits[[disease_id]]$development_metadata),
    held_out_n = nrow(general_splits[[disease_id]]$held_out_metadata),
    development_cases = sum(
      general_splits[[disease_id]]$development_metadata$outcome ==
        POSITIVE_CLASS
    ),
    development_controls = sum(
      general_splits[[disease_id]]$development_metadata$outcome ==
        NEGATIVE_CLASS
    ),
    held_out_cases = sum(
      general_splits[[disease_id]]$held_out_metadata$outcome ==
        POSITIVE_CLASS
    ),
    held_out_controls = sum(
      general_splits[[disease_id]]$held_out_metadata$outcome ==
        NEGATIVE_CLASS
    )
  )

  safe_write_csv(
    split_audit,
    file.path(
      AUDIT_ROOT,
      paste0(
        disease_id,
        "_development_heldout_split_audit.csv"
      )
    )
  )
}


# ==============================================================================
# 21. RUN REQUESTED REPRESENTATION WORKFLOWS
# ==============================================================================

if (run_requested("full")) {
  for (disease_id in names(disease_objects)) {
    run_full_development_representation_stage(
      disease_object = disease_objects[[disease_id]],
      split_object = general_splits[[disease_id]]
    )
  }
}

if (run_requested("experiment1")) {
  run_experiment_1(
    disease_objects = disease_objects
  )
}

if (run_requested("experiment2")) {
  run_experiment_2(
    disease_objects = disease_objects,
    general_splits = general_splits
  )
}

if (run_requested("experiment3")) {
  run_experiment_3(
    disease_objects = disease_objects,
    general_splits = general_splits
  )
}

if (run_requested("experiment4")) {
  run_experiment_4(
    disease_objects = disease_objects,
    general_splits = general_splits
  )
}


# ==============================================================================
# 22. FINAL REPRODUCIBILITY RECORDS
# ==============================================================================

save_config_snapshot(
  CONFIG_PATH,
  SESSION_ROOT
)

write_session_info(
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage01_representation_specification_sessionInfo.txt"
  )
)

run_manifest <- data.table(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  config = CONFIG_PATH,
  modes = paste(RUN_MODES, collapse = ","),
  diseases = paste(REQUESTED_DISEASES, collapse = ","),
  seed = GLOBAL_SEED,
  strict_audits = STRICT_AUDITS,
  overwrite_manifests = OVERWRITE_MANIFESTS
)

safe_write_csv(
  run_manifest,
  file.path(
    SESSION_ROOT,
    "TRACE_Omics_stage01_run_manifest.csv"
  )
)

cat(
  "\n============================================================\n",
  "TRACE-OMICS STAGE 01 COMPLETE\n",
  "Representation specification finished successfully.\n",
  "Output root: ",
  STAGE_ROOT,
  "\n============================================================\n",
  sep = ""
)

# ==============================================================================
# End of scripts/01_representation_specification/01_representation_specification.R
# ==============================================================================

