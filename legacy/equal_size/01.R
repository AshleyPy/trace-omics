# ======================================================================
# TRACE-Omics | Disease-agnostic bulk RNA-seq binary classifier benchmark
# PUBLICATION CANDIDATE (unvalidated; validate against original analyses)
# ======================================================================
# Gene-level raw integer counts (genes x samples) + sample metadata.
# One external YAML config per task: TB, UC, RA, T1D, CD or another case-
# control problem. The disease name is never used as a predictor.
#
# Two sampling modes:
#   equal_size: draw the requested number of cases and controls per task
#   all: use every aligned observation, with class-stratified holdout
#
# Original training-fitted VST/fSVA/ComBat-Ref/ComBat-Seq, five feature
# selectors and seven learners are retained. LOCO requires >= 2 cohorts,
# and transformation feasibility depends on training cohort availability.
#
# WITHOUT EXACT FROZEN PARTITION AND CV FOLD MANIFESTS: a newly sampled
# experiment is NOT a numerical reproduction of a dissertation result.
# Outputs can contain sample-level data and MUST remain PRIVATE.
# ======================================================================

# ======================================================================
# 0. CLEAN ENVIRONMENT
# ======================================================================

# Run in a fresh Rscript session. Do not delete objects from a user's
# interactive/global environment. A separate process also provides a
# more reliable record of the software environment.
gc(verbose = FALSE)


# ======================================================================
# 1. PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "SummarizedExperiment",
  "sva",
  "matrixStats",
  "cluster",
  "neuroCombat",
  "edgeR",
  "mlr3",
  "mlr3filters",
  "mlr3learners",
  "mlr3extralearners",
  "ranger",
  "lightgbm",
  "xgboost",
  "e1071",
  "C50",
  "pROC",
  "ggplot2",
  "RColorBrewer",
  "viridisLite",
  "scales"
)

# Minimal dependencies for data-contract and partition validation.
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install 'data.table' before running validation.")
}
suppressPackageStartupMessages(library(data.table))


# ======================================================================
# 2. SETTINGS
# ======================================================================

SEED <- 7832L
# Record the intended RNG algorithms explicitly (R >= 3.6).
# This complements, but does not replace, version-locked dependencies.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")

POSITIVE_CLASS <- "case"
NEGATIVE_CLASS <- "control"

EQUAL_CASES <- 41L # replaced by YAML (or observed count in all mode)
EQUAL_CONTROLS <- 41L # replaced by YAML (or observed count in all mode)

FINAL_ICE_CASES <- 8L # replaced by YAML / derived holdout fraction
FINAL_ICE_CONTROLS <- 8L # replaced by YAML / derived holdout fraction

TOP_N <- 50L

CV_FOLDS <- 5L
CV_REPEATS <- 1L

MIN_TOTAL_TRAINING_COUNT <- 10L
MAX_SV <- 50L

MAX_PCA_K <- 50L
SILHOUETTE_NEAR_MAX <- 0.95

representation_ids <- c(
  "VST_fSVA",
  "VST_ComBatRef",
  "ComBatSeq_log1p"
)

representation_labels <- c(
  VST_fSVA = "VST + fSVA",
  VST_ComBatRef = "VST + ComBat-Ref",
  ComBatSeq_log1p = "ComBat-Seq + log1p"
)

fsf_ids <- c(
  "information_gain",
  "importance",
  "boruta",
  "mrmr",
  "cmim"
)

fsf_labels <- c(
  information_gain = "Information gain",
  importance = "RF importance",
  boruta = "Boruta",
  mrmr = "mRMR",
  cmim = "CMIM"
)

learner_ids <- c(
  "Naive_Bayes",
  "Linear_SVM",
  "C5_0",
  "Random_Forest",
  "LightGBM",
  "XGBoost",
  "RBF_SVM"
)

learner_labels <- c(
  Naive_Bayes = "Naive Bayes",
  Linear_SVM = "Linear SVM",
  C5_0 = "C5.0",
  Random_Forest = "Random Forest",
  LightGBM = "LightGBM",
  XGBoost = "XGBoost",
  RBF_SVM = "RBF SVM"
)

set.seed(SEED)


# ======================================================================
# 3. SAFE SEED HELPER
# ======================================================================

seed_for <- function(...) {
  values <- as.numeric(c(...))
  weights <- seq_along(values) * 1009
  candidate <- as.numeric(SEED) + sum(values * weights)
  as.integer(candidate %% 2000000000)
}


# ======================================================================
# 4. CONFIGURATION AND PRIVATE INPUT PATHS
# ======================================================================

# Inputs are supplied explicitly; only non-sensitive task parameters go in
# YAML. This avoids hardcoded institution or disease file paths.
usage <- paste(
  "Required: --config FILE.yml --counts FILE --metadata FILE --output-dir NEW_DIR",
  "Optional: --partition-manifest FILE.csv --fold-manifest FILE.csv",
  "Optional: --validate-only true|false (input and design checks only)",
  "Counts: .rds matrix or .csv/.tsv with first column gene_id.",
  "Metadata: .rds data.frame or .csv/.tsv, columns specified in YAML.",
  "Outputs include sensitive sample-level data. Use a private directory OUTSIDE Git.",
  sep = "\n"
)
parse_cli <- function(args) {
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(usage, "\n")
    quit(save = "no", status = 0L)
  }
  if (length(args) == 0L || length(args) %% 2L != 0L) {
    stop(usage, call. = FALSE)
  }
  keys <- args[seq(1L, length(args), by = 2L)]
  values <- args[seq(2L, length(args), by = 2L)]
  allowed <- c("--config", "--counts", "--metadata", "--output-dir",
               "--partition-manifest", "--fold-manifest", "--validate-only")
  if (anyDuplicated(keys) || !all(keys %in% allowed) ||
      any(!nzchar(values)) ||
      !all(c("--config", "--counts", "--metadata", "--output-dir") %in% keys)) {
    stop("Invalid CLI arguments.\n", usage, call. = FALSE)
  }
  stats::setNames(as.list(values), sub("^--", "", keys))
}
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install R package 'yaml' to read task configuration.")
}
config_path <- normalizePath(cli[["config"]], winslash = "/", mustWork = TRUE)
config <- yaml::read_yaml(config_path)
if (!is.list(config) || !is.list(config$design) || !is.list(config$columns)) {
  stop("YAML must define 'dataset_id', 'columns' and 'design'.")
}
DATASET_ID <- as.character(config$dataset_id)
if (length(DATASET_ID) != 1L || is.na(DATASET_ID) ||
    !grepl("^[A-Za-z][A-Za-z0-9_-]{0,40}$", DATASET_ID)) {
  stop("dataset_id must be a short alphanumeric label (e.g. TB, UC, RA, T1D).")
}
SAMPLING_MODE <- as.character(config$design$sampling_mode)
if (!SAMPLING_MODE %in% c("equal_size", "all")) {
  stop("design.sampling_mode must be equal_size or all.")
}
get_int <- function(x, name, min_value = 1L) {
  n <- suppressWarnings(as.numeric(x))
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n != floor(n) ||
      n < min_value || n > .Machine$integer.max) {
    stop("Invalid integer setting: ", name)
  }
  as.integer(n)
}
SEED <- get_int(config$design$seed, "design.seed", 0L)
CV_FOLDS <- get_int(config$design$cv_folds, "design.cv_folds", 2L)
CV_REPEATS <- get_int(config$design$cv_repeats, "design.cv_repeats")
TOP_N <- get_int(config$design$top_genes, "design.top_genes")
MAX_SV <- get_int(config$design$max_surrogate_variables,
                  "design.max_surrogate_variables", 0L)
if (SAMPLING_MODE == "equal_size") {
  EQUAL_CASES <- get_int(config$design$equal_cases, "design.equal_cases")
  EQUAL_CONTROLS <- get_int(config$design$equal_controls, "design.equal_controls")
  FINAL_ICE_CASES <- get_int(config$design$holdout_cases, "design.holdout_cases")
  FINAL_ICE_CONTROLS <- get_int(config$design$holdout_controls,
                                "design.holdout_controls")
  if (FINAL_ICE_CASES >= EQUAL_CASES ||
      FINAL_ICE_CONTROLS >= EQUAL_CONTROLS) {
    stop("Holdout class counts must be smaller than selected class counts.")
  }
} else {
  fraction <- as.numeric(config$design$holdout_fraction)
  if (length(fraction) != 1L || is.na(fraction) ||
      fraction <= 0 || fraction >= 0.5) {
    stop("For all mode, holdout_fraction must be in (0, 0.5).")
  }
}
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
set.seed(SEED)

# Mappings are required; never infer the target from arbitrary column names.
cols <- config$columns
needed_cols <- c("sample_id", "cohort", "outcome")
if (!all(needed_cols %in% names(cols))) {
  stop("columns must specify sample_id, cohort and outcome.")
}
for (field in needed_cols) {
  if (!is.character(cols[[field]]) || length(cols[[field]]) != 1L ||
      is.na(cols[[field]]) || !nzchar(cols[[field]])) {
    stop("Invalid column name for: ", field)
  }
}
if (anyDuplicated(unlist(cols[needed_cols], use.names = FALSE))) {
  stop("Configured sample_id, cohort and outcome source columns must differ.")
}
if (!is.null(cols$participant_id) &&
    (!is.character(cols$participant_id) || length(cols$participant_id) != 1L)) {
  stop("Invalid columns.participant_id")
}
CASE_LABELS <- unique(tolower(trimws(c("case", as.character(config$labels$case)))))
CONTROL_LABELS <- unique(tolower(trimws(c("control", as.character(config$labels$control)))))
if (anyNA(CASE_LABELS) || anyNA(CONTROL_LABELS) ||
    any(CASE_LABELS == "") || any(CONTROL_LABELS == "") ||
    length(intersect(CASE_LABELS, CONTROL_LABELS))) {
  stop("Invalid or overlapping case/control label aliases.")
}
`%||%` <- function(x, default) if (is.null(x)) default else x
VALIDATE_ONLY <- identical(tolower(cli[["validate-only"]] %||% "false"), "true")

counts_path <- normalizePath(cli[["counts"]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(cli[["metadata"]], winslash = "/", mustWork = TRUE)
partition_manifest_path <- cli[["partition-manifest"]]
fold_manifest_path <- cli[["fold-manifest"]]
if (!is.null(partition_manifest_path)) {
  partition_manifest_path <- normalizePath(partition_manifest_path,
                                           winslash = "/", mustWork = TRUE)
}
if (!is.null(fold_manifest_path)) {
  fold_manifest_path <- normalizePath(fold_manifest_path,
                                      winslash = "/", mustWork = TRUE)
}
if (xor(is.null(partition_manifest_path), is.null(fold_manifest_path))) {
  warning("Only one frozen manifest supplied. Exact partition/fold reproduction ",
          "requires BOTH manifests.", call. = FALSE)
}
output_arg <- path.expand(cli[["output-dir"]])
output_parent <- normalizePath(dirname(output_arg),
                               winslash = "/", mustWork = TRUE)
output_root <- file.path(output_parent, basename(output_arg))
inside_git_tree <- function(start_dir) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, ".git")) ||
        dir.exists(file.path(current, ".git"))) return(TRUE)
    parent <- dirname(current)
    if (identical(parent, current)) return(FALSE)
    current <- parent
  }
}
if (inside_git_tree(output_parent)) {
  stop("Refusing to write research outputs inside a Git repository. ",
       "Choose an external private directory.")
}
if (file.exists(output_root) || dir.exists(output_root)) {
  stop("Output path already exists; specify a new directory: ", output_root)
}
analysis_mode <- if (!is.null(partition_manifest_path) &&
                     !is.null(fold_manifest_path)) {
  "frozen_manifest_reproduction_candidate"
} else {
  "new_seeded_replication_or_partial_manifest"
}
message("Dataset: ", DATASET_ID, " | sampling: ", SAMPLING_MODE,
        " | analysis mode: ", analysis_mode)
message("All outputs are RESTRICTED and must remain private: ", output_root)

audit_dir <- file.path(
  output_root,
  "00_audits"
)

split_dir <- file.path(
  output_root,
  "01_split"
)

diagnostic_dir <- file.path(
  output_root,
  "02_diagnostics"
)

concordance_dir <- file.path(
  output_root,
  "03_feature_concordance"
)

benchmark_dir <- file.path(
  output_root,
  "04_benchmark"
)

ccv_dir <- file.path(
  output_root,
  "05_CCV"
)

final_ice_dir <- file.path(
  output_root,
  "99_FINAL_ICE"
)

for (current_dir in c(
  output_root,
  audit_dir,
  split_dir,
  diagnostic_dir,
  concordance_dir,
  benchmark_dir,
  ccv_dir,
  final_ice_dir
)) {
  dir.create(
    current_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ======================================================================
# 5. STANDARDISE BINARY OUTCOME
# ======================================================================

standardise_outcome <- function(outcome) {
  values <- tolower(trimws(as.character(outcome)))
  values[values %in% CASE_LABELS] <- "case"
  values[values %in% CONTROL_LABELS] <- "control"
  invalid <- unique(as.character(outcome)[!values %in% c("case", "control") |
                                           is.na(values)])
  if (length(invalid)) {
    stop("Unrecognised outcome labels: ", paste(invalid, collapse = ", "),
         ". Add explicit aliases in config.labels.")
  }
  factor(values, levels = c("control", "case"))
}


# ======================================================================
# 6. LOAD + ALIGN GENE-LEVEL RAW DATA
# ======================================================================

read_count_matrix <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") return(as.matrix(readRDS(path)))
  if (!extension %in% c("csv", "tsv")) {
    stop("Counts must be .rds, .csv or .tsv with genes in rows.")
  }
  sep <- if (extension == "csv") "," else "\t"
  d <- data.table::fread(path, sep = sep, check.names = FALSE)
  if (ncol(d) < 3L) stop("Counts table needs gene_id plus >= 2 samples.")
  genes <- as.character(d[[1L]])
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("First counts column must contain unique nonblank gene IDs.")
  }
  d[, (names(d)[1L]) := NULL]
  m <- as.matrix(d)
  if (!is.numeric(m)) stop("Count columns must be numeric; avoid quoted strings.")
  rownames(m) <- genes
  m
}
read_sample_metadata <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") return(data.table::as.data.table(readRDS(path)))
  if (extension %in% c("csv", "tsv")) {
    return(data.table::fread(path, sep = if (extension == "csv") "," else "\t"))
  }
  stop("Metadata must be .rds, .csv or .tsv.")
}
rna_counts <- read_count_matrix(counts_path)
rna_metadata <- read_sample_metadata(metadata_path)

if (
  is.null(rownames(rna_counts)) ||
    is.null(colnames(rna_counts))
) {
  stop(
    "Raw-count matrix requires gene and sample identifiers."
  )
}

storage.mode(rna_counts) <- "numeric"

if (
  anyNA(rna_counts) ||
    !all(is.finite(rna_counts)) ||
    any(rna_counts < 0)
) {
  stop(
    "Raw counts contain invalid values."
  )
}

if (any(
  abs(
    rna_counts -
      round(rna_counts)
  ) > 1e-8
)) {
  stop(
    "RNA-seq counts must be non-negative integers (not VST, TPM or CPM)."
  )
}

if (any(rna_counts > .Machine$integer.max)) {
  stop("Counts exceed R's 32-bit integer range; check input encoding.")
}
rna_counts <- round(rna_counts)
storage.mode(rna_counts) <- "integer"

if (anyDuplicated(rownames(rna_counts))) {
  stop("Duplicated gene identifiers detected.")
}

if (anyDuplicated(colnames(rna_counts))) {
  stop("Duplicated count-matrix sample identifiers detected.")
}

# Map task-specific source fields to the framework's internal data contract.
# These provenance fields are never included among predictive features.
for (field in needed_cols) {
  if (!cols[[field]] %in% names(rna_metadata)) {
    stop("Metadata missing configured column: ", cols[[field]])
  }
}
if (!is.null(cols$participant_id) &&
    !cols$participant_id %in% names(rna_metadata)) {
  stop("Metadata missing configured participant column: ", cols$participant_id)
}
rna_metadata[, sample_id := trimws(as.character(get(cols$sample_id)))]
rna_metadata[, cohort_group := trimws(as.character(get(cols$cohort)))]
rna_metadata[, outcome := standardise_outcome(get(cols$outcome))]
if (anyNA(rna_metadata$sample_id) || any(!nzchar(rna_metadata$sample_id))) {
  stop("Sample identifiers contain blanks or missing values.")
}
if (anyNA(rna_metadata$cohort_group) || any(!nzchar(rna_metadata$cohort_group))) {
  stop("Cohort identifiers contain blanks or missing values.")
}
if (!is.null(cols$participant_id)) {
  rna_metadata[, participant_id := trimws(as.character(get(cols$participant_id)))]
}
if (anyDuplicated(rna_metadata$sample_id)) {
  stop("Duplicated metadata sample IDs detected.")
}

common_samples <- colnames(rna_counts)[
  colnames(rna_counts) %in%
    rna_metadata$sample_id
]

if (length(common_samples) < 4L) {
  stop("At least four aligned observations are required.")
}

rna_counts <- rna_counts[
  ,
  common_samples,
  drop = FALSE
]

metadata_match <- match(
  common_samples,
  rna_metadata$sample_id
)

if (anyNA(metadata_match)) {
  stop(
    "Count/metadata alignment failed."
  )
}

rna_metadata <- rna_metadata[
  metadata_match
]

stopifnot(
  identical(
    colnames(rna_counts),
    rna_metadata$sample_id
  )
)

# Counts for all mode are taken from the ALIGNED input population.
# Its natural class imbalance is retained rather than forced to 50:50.
if (SAMPLING_MODE == "all") {
  EQUAL_CASES <- sum(rna_metadata$outcome == POSITIVE_CLASS)
  EQUAL_CONTROLS <- sum(rna_metadata$outcome == NEGATIVE_CLASS)
  FINAL_ICE_CASES <- max(1L, as.integer(round(EQUAL_CASES * fraction)))
  FINAL_ICE_CONTROLS <- max(1L, as.integer(round(EQUAL_CONTROLS * fraction)))
}
if (EQUAL_CASES <= FINAL_ICE_CASES ||
    EQUAL_CONTROLS <= FINAL_ICE_CONTROLS ||
    sum(rna_metadata$outcome == POSITIVE_CLASS) < EQUAL_CASES ||
    sum(rna_metadata$outcome == NEGATIVE_CLASS) < EQUAL_CONTROLS) {
  stop("Insufficient cases/controls for the configured selection and holdout.")
}
if (length(unique(rna_metadata$cohort_group)) < 2L) {
  stop("This LOCO-enabled benchmark requires >= 2 cohorts. ",
       "A single-cohort task needs a separate internal-validation workflow.")
}
rna_metadata[, disease := DATASET_ID]


# Optional participant-level leakage audits are enforced if participant IDs
# are supplied. Study-scoped keys prevent collisions across source studies.
if ("participant_id" %in% names(rna_metadata)) {
  rna_metadata[, participant_id := trimws(as.character(participant_id))]
  if (anyNA(rna_metadata$participant_id) ||
      any(!nzchar(rna_metadata$participant_id))) {
    stop("participant_id is incomplete; cannot verify participant independence.")
  }
  rna_metadata[, participant_key := paste(cohort_group, participant_id,
                                         sep = "::")]
} else {
  warning("participant_id absent: participant-level split independence ",
          "cannot be independently verified.", call. = FALSE)
}

# File fingerprints identify the exact external inputs; MD5 is used solely
# as an integrity fingerprint, not as a security or anonymisation measure.
input_fingerprints <- data.frame(
  role = c("configuration", "counts", "metadata", "partition_manifest", "fold_manifest"),
  path = c(config_path, counts_path, metadata_path,
           if (is.null(partition_manifest_path)) NA_character_ else partition_manifest_path,
           if (is.null(fold_manifest_path)) NA_character_ else fold_manifest_path),
  stringsAsFactors = FALSE
)
input_fingerprints$md5 <- vapply(
  input_fingerprints$path,
  function(p) if (is.na(p)) NA_character_ else unname(tools::md5sum(p)),
  character(1)
)
utils::write.csv(input_fingerprints,
                 file.path(audit_dir, paste0(DATASET_ID, "_private_input_fingerprints.csv")),
                 row.names = FALSE)


observed_cohorts <- sort(
  unique(
    as.character(
      rna_metadata$cohort_group
    )
  )
)

cat(
  paste0("\nObserved ", DATASET_ID, " cohorts after alignment:\n"),
  paste(
    observed_cohorts,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


raw_audit <- data.table::data.table(
  Population = paste0("Full ", DATASET_ID),
  n = ncol(rna_counts),
  Controls = sum(
    rna_metadata$outcome ==
      NEGATIVE_CLASS
  ),
  Cases = sum(
    rna_metadata$outcome ==
      POSITIVE_CLASS
  ),
  Cohorts = data.table::uniqueN(
    rna_metadata$cohort_group
  ),
  Genes = nrow(rna_counts)
)

data.table::fwrite(
  raw_audit,
  file.path(
    audit_dir,
    paste0(DATASET_ID, "_full_raw_input_audit.csv")
  )
)

cat(
  "\n============================================================\n",
  paste0("FULL ", DATASET_ID, " RAW DATA\n"),
  "============================================================\n",
  sep = ""
)

print(
  raw_audit
)

cat(
  paste0("\nFull ", DATASET_ID, " cohort x outcome:\n")
)

print(
  table(
    rna_metadata$cohort_group,
    rna_metadata$outcome
  )
)


# ======================================================================
# 7. PROPORTIONAL ALLOCATION HELPER
# ======================================================================

allocate_integer_target <- function(
  counts_by_stratum,
  target,
  retain_at_least_one = FALSE,
  give_each_if_possible = TRUE
) {

  stratum_names <- names(
    counts_by_stratum
  )

  counts_by_stratum <- as.integer(
    counts_by_stratum
  )

  names(counts_by_stratum) <- stratum_names

  if (target < 0L) {
    stop("Allocation target must be non-negative.")
  }

  if (sum(counts_by_stratum) < target) {
    stop(
      "Allocation target exceeds available observations."
    )
  }

  eligible <- counts_by_stratum > 0L

  allocation <- integer(
    length(counts_by_stratum)
  )

  names(allocation) <- names(
    counts_by_stratum
  )

  capacity <- counts_by_stratum

  if (retain_at_least_one) {
    capacity <- pmax(
      counts_by_stratum -
        as.integer(eligible),
      0L
    )
  }

  if (sum(capacity) < target) {
    stop(
      "Cannot allocate requested target while retaining at least one observation in each non-empty stratum."
    )
  }

  if (
    give_each_if_possible &&
      target >=
        sum(
          eligible &
            capacity > 0L
        )
  ) {

    initial <- which(
      eligible &
        capacity > 0L
    )

    allocation[
      initial
    ] <- 1L
  }

  allocation <- pmin(
    allocation,
    capacity
  )

  remaining <- target -
    sum(allocation)

  while (remaining > 0L) {

    available <- which(
      allocation <
        capacity
    )

    if (length(available) == 0L) {
      stop(
        "Allocation failed before reaching target."
      )
    }

    residual_capacity <- capacity[
      available
    ] -
      allocation[
        available
      ]

    proportional_weight <- counts_by_stratum[
      available
    ]

    if (sum(proportional_weight) <= 0) {
      proportional_weight <- residual_capacity
    }

    raw_add <- remaining *
      proportional_weight /
      sum(proportional_weight)

    floor_add <- pmin(
      floor(raw_add),
      residual_capacity
    )

    if (sum(floor_add) > 0L) {

      allocation[
        available
      ] <- allocation[
        available
      ] +
        floor_add

      remaining <- target -
        sum(allocation)

      if (remaining <= 0L) {
        break
      }
    }

    available <- which(
      allocation <
        capacity
    )

    if (length(available) == 0L) {
      break
    }

    fractional <- numeric(
      length(available)
    )

    for (index in seq_along(available)) {

      stratum_index <- available[[index]]

      fractional[[index]] <-
        (
          target *
            counts_by_stratum[[stratum_index]] /
            sum(counts_by_stratum)
        ) -
        floor(
          target *
            counts_by_stratum[[stratum_index]] /
            sum(counts_by_stratum)
        )
    }

    order_add <- order(
      fractional,
      counts_by_stratum[
        available
      ],
      decreasing = TRUE
    )

    for (ordered_index in order_add) {

      if (remaining <= 0L) {
        break
      }

      stratum_index <- available[
        ordered_index
      ]

      if (
        allocation[[stratum_index]] <
          capacity[[stratum_index]]
      ) {

        allocation[[stratum_index]] <-
          allocation[[stratum_index]] +
          1L

        remaining <- remaining -
          1L
      }
    }
  }

  if (sum(allocation) != target) {
    stop(
      "Final allocation does not equal requested target."
    )
  }

  allocation
}


# ======================================================================
# 8. CONFIGURABLE COHORT-STRATIFIED SAMPLE SELECTION
# ======================================================================

select_equal_size_ids <- function(
  metadata,
  outcome_value,
  target,
  seed
) {

  current <- metadata[
    outcome ==
      outcome_value
  ]

  cohort_counts <- table(
    current$cohort_group
  )

  allocation <- allocate_integer_target(
    counts_by_stratum =
      cohort_counts,
    target =
      target,
    retain_at_least_one =
      FALSE,
    give_each_if_possible =
      TRUE
  )

  selected_ids <- character(0)

  allocation_rows <- list()

  row_counter <- 0L

  for (cohort_name in names(allocation)) {

    requested <- allocation[[cohort_name]]

    cohort_ids <- current[
      cohort_group ==
        cohort_name,
      sample_id
    ]

    set.seed(
      seed_for(
        seed,
        match(
          cohort_name,
          names(allocation)
        ),
        ifelse(
          outcome_value ==
            POSITIVE_CLASS,
          1,
          2
        )
      )
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

    selected_ids <- c(
      selected_ids,
      chosen
    )

    row_counter <- row_counter + 1L

    allocation_rows[[row_counter]] <-
      data.table::data.table(
        Outcome = outcome_value,
        Cohort = cohort_name,
        Available = length(cohort_ids),
        Selected = requested
      )
  }

  list(
    ids = selected_ids,
    allocation = data.table::rbindlist(
      allocation_rows
    )
  )
}

# A frozen manifest prevents an ostensibly identical rerun from selecting
# different observations because of input order, R version or RNG behaviour.
# It is PRIVATE because it records sample identifiers and case/control status.
frozen_partition <- NULL
if (!is.null(partition_manifest_path)) {
  frozen_partition <- data.table::fread(partition_manifest_path)
  needed <- c("sample_id", "partition")
  if (!all(needed %in% names(frozen_partition))) {
    stop("Partition manifest must contain: ", paste(needed, collapse = ", "))
  }
  frozen_partition[, sample_id := as.character(sample_id)]
  frozen_partition[, partition := as.character(partition)]
  if (anyNA(frozen_partition$sample_id) ||
      any(!nzchar(frozen_partition$sample_id)) ||
      anyDuplicated(frozen_partition$sample_id) ||
      !all(frozen_partition$partition %in% c("development", "final_ice")) ||
      nrow(frozen_partition) != EQUAL_CASES + EQUAL_CONTROLS ||
      !all(frozen_partition$sample_id %in% rna_metadata$sample_id)) {
    stop("Partition manifest has invalid, missing, duplicated or ",
         "unrecognised sample IDs/partition labels.")
  }
  if ("selection_order" %in% names(frozen_partition)) {
    if (anyNA(frozen_partition$selection_order) ||
        anyDuplicated(frozen_partition$selection_order)) {
      stop("Invalid selection_order in partition manifest.")
    }
    data.table::setorder(frozen_partition, selection_order)
  }
  matched <- rna_metadata[match(frozen_partition$sample_id, sample_id)]
  if ("outcome" %in% names(frozen_partition) &&
      !identical(as.character(frozen_partition$outcome),
                 as.character(matched$outcome))) {
    stop("Partition manifest outcome labels disagree with input metadata.")
  }
  if ("cohort_group" %in% names(frozen_partition) &&
      !identical(as.character(frozen_partition$cohort_group),
                 as.character(matched$cohort_group))) {
    stop("Partition manifest cohort labels disagree with input metadata.")
  }
  if (sum(matched$outcome == POSITIVE_CLASS) != EQUAL_CASES ||
      sum(matched$outcome == NEGATIVE_CLASS) != EQUAL_CONTROLS ||
      sum(frozen_partition$partition == "final_ice" &
          matched$outcome == POSITIVE_CLASS) != FINAL_ICE_CASES ||
      sum(frozen_partition$partition == "final_ice" &
          matched$outcome == NEGATIVE_CLASS) != FINAL_ICE_CONTROLS) {
    stop("Frozen partition does not match configured class and holdout sizes.")
  }
}

manifest_selection <- function(label) {
  available <- rna_metadata[outcome == label, .N, by = cohort_group]
  data.table::setnames(available, "N", "Available")
  chosen <- rna_metadata[
    sample_id %in% frozen_partition$sample_id & outcome == label,
    .N, by = cohort_group]
  data.table::setnames(chosen, "N", "Selected")
  allocation <- merge(available, chosen, by = "cohort_group", all.x = TRUE)
  allocation[is.na(Selected), Selected := 0L]
  data.table::setnames(allocation, "cohort_group", "Cohort")
  allocation[, Outcome := label]
  data.table::setcolorder(allocation, c("Outcome", "Cohort", "Available", "Selected"))
  list(ids = frozen_partition$sample_id[
         frozen_partition$sample_id %in%
           rna_metadata[outcome == label, sample_id]],
       allocation = allocation)
}

if (is.null(frozen_partition)) {
  selected_cases <- select_equal_size_ids(
    metadata = rna_metadata,
    outcome_value = POSITIVE_CLASS,
    target = EQUAL_CASES,
    seed = SEED + 10
  )
  selected_controls <- select_equal_size_ids(
    metadata = rna_metadata,
    outcome_value = NEGATIVE_CLASS,
    target = EQUAL_CONTROLS,
    seed = SEED + 20
  )
  equal_ids <- c(selected_cases$ids, selected_controls$ids)
} else {
  selected_cases <- manifest_selection(POSITIVE_CLASS)
  selected_controls <- manifest_selection(NEGATIVE_CLASS)
  equal_ids <- frozen_partition$sample_id
}

equal_size_audit <- data.table::data.table(
  Requested_cases = EQUAL_CASES,
  Selected_cases = length(selected_cases$ids),
  Requested_controls = EQUAL_CONTROLS,
  Selected_controls = length(selected_controls$ids),
  Total_selected = length(equal_ids),
  Duplicated_IDs = anyDuplicated(equal_ids)
)

print(
  equal_size_audit
)

if (
  length(selected_cases$ids) != EQUAL_CASES ||
    length(selected_controls$ids) != EQUAL_CONTROLS ||
    length(equal_ids) !=
      EQUAL_CASES +
        EQUAL_CONTROLS ||
    anyDuplicated(equal_ids)
) {
  stop(
    paste0(
      "Sample selection failed. ",
      "Cases selected = ",
      length(selected_cases$ids),
      "/",
      EQUAL_CASES,
      "; controls selected = ",
      length(selected_controls$ids),
      "/",
      EQUAL_CONTROLS,
      "; total = ",
      length(equal_ids),
      "; duplicated IDs = ",
      anyDuplicated(equal_ids),
      ". Inspect the dataset-specific cohort allocation file."
    )
  )
}

equal_metadata <- rna_metadata[
  match(
    equal_ids,
    sample_id
  )
]

equal_counts <- rna_counts[
  ,
  equal_metadata$sample_id,
  drop = FALSE
]

stopifnot(
  sum(
    equal_metadata$outcome ==
      POSITIVE_CLASS
  ) == EQUAL_CASES,
  sum(
    equal_metadata$outcome ==
      NEGATIVE_CLASS
  ) == EQUAL_CONTROLS
)

equal_allocation <- data.table::rbindlist(
  list(
    selected_cases$allocation,
    selected_controls$allocation
  )
)

data.table::fwrite(
  equal_allocation,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_selected_cohort_allocation.csv")
  )
)


# ======================================================================
# 9. FREEZE CONFIGURED DEVELOPMENT / HOLDOUT SPLIT
# ======================================================================

select_final_ice_ids <- function(
  metadata,
  outcome_value,
  target,
  seed
) {

  current <- metadata[
    outcome ==
      outcome_value
  ]

  cohort_counts <- table(
    current$cohort_group
  )

  allocation <- allocate_integer_target(
    counts_by_stratum =
      cohort_counts,
    target =
      target,
    retain_at_least_one =
      TRUE,
    give_each_if_possible =
      TRUE
  )

  test_ids <- character(0)

  rows <- list()
  row_counter <- 0L

  for (cohort_name in names(allocation)) {

    requested <- allocation[[cohort_name]]

    cohort_ids <- current[
      cohort_group ==
        cohort_name,
      sample_id
    ]

    set.seed(
      seed_for(
        seed,
        match(
          cohort_name,
          names(allocation)
        ),
        ifelse(
          outcome_value ==
            POSITIVE_CLASS,
          11,
          22
        )
      )
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

    test_ids <- c(
      test_ids,
      chosen
    )

    row_counter <- row_counter + 1L

    rows[[row_counter]] <-
      data.table::data.table(
        Outcome = outcome_value,
        Cohort = cohort_name,
        Equal_size_available = length(cohort_ids),
        Final_ice_selected = requested,
        Development_retained =
          length(cohort_ids) -
          requested
      )
  }

  list(
    ids = test_ids,
    allocation = data.table::rbindlist(
      rows
    )
  )
}

if (is.null(frozen_partition)) {
  ice_cases <- select_final_ice_ids(
    metadata = equal_metadata,
    outcome_value = POSITIVE_CLASS,
    target = FINAL_ICE_CASES,
    seed = SEED + 30
  )
  ice_controls <- select_final_ice_ids(
    metadata = equal_metadata,
    outcome_value = NEGATIVE_CLASS,
    target = FINAL_ICE_CONTROLS,
    seed = SEED + 40
  )
  final_ice_ids <- c(ice_cases$ids, ice_controls$ids)
} else {
  final_ice_ids <- frozen_partition[partition == "final_ice", sample_id]
  # These IDs were validated against the configured holdout design.
}

# Persist the ACTUAL selected IDs and partitions. Keep this file private;
# it is essential for exact numerical reproduction of the sampling design.
private_partition_manifest <- data.table::data.table(
  sample_id = equal_metadata$sample_id,
  cohort_group = equal_metadata$cohort_group,
  outcome = as.character(equal_metadata$outcome),
  partition = ifelse(equal_metadata$sample_id %in% final_ice_ids,
                     "final_ice", "development"),
  selection_order = seq_len(nrow(equal_metadata))
)
data.table::fwrite(private_partition_manifest,
  file.path(split_dir, paste0(DATASET_ID, "_frozen_partition_manifest_PRIVATE.csv")))


development_ids <- setdiff(
  equal_metadata$sample_id,
  final_ice_ids
)

development_metadata <- equal_metadata[
  match(
    development_ids,
    sample_id
  )
]

final_ice_metadata <- equal_metadata[
  match(
    final_ice_ids,
    sample_id
  )
]

development_counts <- equal_counts[
  ,
  development_metadata$sample_id,
  drop = FALSE
]

final_ice_counts <- equal_counts[
  ,
  final_ice_metadata$sample_id,
  drop = FALSE
]

if ("participant_key" %in% names(equal_metadata)) {
  shared_participants <- intersect(development_metadata$participant_key,
                                   final_ice_metadata$participant_key)
  if (length(shared_participants) > 0L) {
    stop("Participant leakage: participant appears in both development ",
         "and final holdout. Supply participant-grouped frozen partitions.")
  }
}

stopifnot(
  ncol(development_counts) == EQUAL_CASES + EQUAL_CONTROLS - FINAL_ICE_CASES - FINAL_ICE_CONTROLS,
  ncol(final_ice_counts) == FINAL_ICE_CASES + FINAL_ICE_CONTROLS,
  sum(
    development_metadata$outcome ==
      POSITIVE_CLASS
  ) == EQUAL_CASES - FINAL_ICE_CASES,
  sum(
    development_metadata$outcome ==
      NEGATIVE_CLASS
  ) == EQUAL_CONTROLS - FINAL_ICE_CONTROLS,
  sum(
    final_ice_metadata$outcome ==
      POSITIVE_CLASS
  ) == FINAL_ICE_CASES,
  sum(
    final_ice_metadata$outcome ==
      NEGATIVE_CLASS
  ) == FINAL_ICE_CONTROLS,
  length(
    intersect(
      development_metadata$sample_id,
      final_ice_metadata$sample_id
    )
  ) == 0L,
  all(
    unique(
      final_ice_metadata$cohort_group
    ) %in%
      unique(
        development_metadata$cohort_group
      )
  )
)

# Strict held-out cohort AUROC requires both classes in every development
# cohort. Fail early rather than silently report undefined metrics.
cohort_class_counts <- table(development_metadata$cohort_group,
                             development_metadata$outcome)
if (any(cohort_class_counts == 0L)) {
  stop("Every development cohort must contain cases and controls for strict LOCO AUROC. ",
       "Adjust the sample allocation or input design.")
}
if (any(table(development_metadata$outcome) < CV_FOLDS)) {
  stop("Each development class requires at least cv_folds observations.")
}

split_summary <- data.table::rbindlist(
  list(
    data.table::data.table(
      Partition = "Selected population",
      n = nrow(equal_metadata),
      Controls = sum(
        equal_metadata$outcome ==
          NEGATIVE_CLASS
      ),
      Cases = sum(
        equal_metadata$outcome ==
          POSITIVE_CLASS
      ),
      Cohorts = data.table::uniqueN(
        equal_metadata$cohort_group
      ),
      IR = max(sum(equal_metadata$outcome == POSITIVE_CLASS), sum(equal_metadata$outcome == NEGATIVE_CLASS)) / min(sum(equal_metadata$outcome == POSITIVE_CLASS), sum(equal_metadata$outcome == NEGATIVE_CLASS))
    ),
    data.table::data.table(
      Partition = "Development",
      n = nrow(development_metadata),
      Controls = sum(
        development_metadata$outcome ==
          NEGATIVE_CLASS
      ),
      Cases = sum(
        development_metadata$outcome ==
          POSITIVE_CLASS
      ),
      Cohorts = data.table::uniqueN(
        development_metadata$cohort_group
      ),
      IR = max(sum(development_metadata$outcome == POSITIVE_CLASS), sum(development_metadata$outcome == NEGATIVE_CLASS)) / min(sum(development_metadata$outcome == POSITIVE_CLASS), sum(development_metadata$outcome == NEGATIVE_CLASS))
    ),
    data.table::data.table(
      Partition = "Final holdout",
      n = nrow(final_ice_metadata),
      Controls = sum(
        final_ice_metadata$outcome ==
          NEGATIVE_CLASS
      ),
      Cases = sum(
        final_ice_metadata$outcome ==
          POSITIVE_CLASS
      ),
      Cohorts = data.table::uniqueN(
        final_ice_metadata$cohort_group
      ),
      IR = max(sum(final_ice_metadata$outcome == POSITIVE_CLASS), sum(final_ice_metadata$outcome == NEGATIVE_CLASS)) / min(sum(final_ice_metadata$outcome == POSITIVE_CLASS), sum(final_ice_metadata$outcome == NEGATIVE_CLASS))
    )
  )
)

cohort_names <- sort(
  unique(
    equal_metadata$cohort_group
  )
)

split_by_cohort_rows <- list()

for (cohort_name in cohort_names) {

  selected_current <- equal_metadata[
    cohort_group ==
      cohort_name
  ]

  train_current <- development_metadata[
    cohort_group ==
      cohort_name
  ]

  ice_current <- final_ice_metadata[
    cohort_group ==
      cohort_name
  ]

  split_by_cohort_rows[[cohort_name]] <-
    data.table::data.table(
      Cohort = cohort_name,
      Selected_n = nrow(selected_current),
      Selected_controls = sum(
        selected_current$outcome ==
          NEGATIVE_CLASS
      ),
      Selected_cases = sum(
        selected_current$outcome ==
          POSITIVE_CLASS
      ),
      Train_n = nrow(train_current),
      Train_controls = sum(
        train_current$outcome ==
          NEGATIVE_CLASS
      ),
      Train_cases = sum(
        train_current$outcome ==
          POSITIVE_CLASS
      ),
      Final_ice_n = nrow(ice_current),
      Final_ice_controls = sum(
        ice_current$outcome ==
          NEGATIVE_CLASS
      ),
      Final_ice_cases = sum(
        ice_current$outcome ==
          POSITIVE_CLASS
      )
    )
}

split_by_cohort <- data.table::rbindlist(
  split_by_cohort_rows
)

data.table::fwrite(
  split_summary,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_train_test_split_summary.csv")
  )
)

data.table::fwrite(
  split_by_cohort,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_train_test_split_by_cohort.csv")
  )
)

data.table::fwrite(
  development_metadata,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_development_metadata.csv")
  )
)

data.table::fwrite(
  final_ice_metadata,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_FINAL_ICE_metadata.csv")
  )
)

saveRDS(
  development_counts,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_development_raw_counts.rds")
  )
)

saveRDS(
  final_ice_counts,
  file.path(
    split_dir,
    paste0(DATASET_ID, "_FINAL_ICE_raw_counts.rds")
  )
)

cat(
  "\n============================================================\n",
  "TRAIN / FINAL-ICE SPLIT\n",
  "============================================================\n",
  sep = ""
)

print(
  split_summary
)

cat(
  "\nBy cohort:\n"
)

print(
  split_by_cohort
)


# ======================================================================
# 10. FOLD MAP: COHORT x OUTCOME STRATIFIED
# ======================================================================

make_stratified_fold_map <- function(
  metadata,
  folds,
  repeats,
  seed
) {

  output_rows <- list()
  counter <- 0L

  for (repeat_index in seq_len(repeats)) {

    strata <- interaction(
      metadata$cohort_group,
      metadata$outcome,
      drop = TRUE,
      lex.order = TRUE
    )

    fold_assignment <- integer(
      nrow(metadata)
    )

    for (stratum_name in levels(strata)) {

      current_index <- which(
        strata ==
          stratum_name
      )

      set.seed(
        seed_for(
          seed,
          repeat_index,
          match(
            stratum_name,
            levels(strata)
          )
        )
      )

      current_index <- sample(
        current_index,
        length(current_index),
        replace = FALSE
      )

      fold_values <- rep(
        seq_len(folds),
        length.out = length(current_index)
      )

      fold_assignment[
        current_index
      ] <- fold_values
    }

    for (sample_index in seq_len(nrow(metadata))) {

      counter <- counter + 1L

      output_rows[[counter]] <-
        data.table::data.table(
          sample_id = metadata$sample_id[[sample_index]],
          cohort_group = metadata$cohort_group[[sample_index]],
          outcome = as.character(
            metadata$outcome[[sample_index]]
          ),
          Repeat = repeat_index,
          Fold = fold_assignment[[sample_index]]
        )
    }
  }

  output <- data.table::rbindlist(
    output_rows
  )

  for (repeat_index in seq_len(repeats)) {

    for (fold_index in seq_len(folds)) {

      validation <- output[
        Repeat ==
          repeat_index &
          Fold ==
            fold_index
      ]

      train_ids <- setdiff(
        metadata$sample_id,
        validation$sample_id
      )

      train_metadata <- metadata[
        sample_id %in%
          train_ids
      ]

      if ("participant_key" %in% names(metadata)) {
        validation_participants <- metadata[
          sample_id %in% validation$sample_id, participant_key]
        if (length(intersect(train_metadata$participant_key,
                             validation_participants)) > 0L) {
          stop("Participant leakage across development CV folds; ",
               "use a participant-grouped fold manifest.")
        }
      }

      if (data.table::uniqueN(
        train_metadata$outcome
      ) < 2L) {
        stop(
          "CV fold has only one training outcome class."
        )
      }

      if (!all(
        unique(
          validation$cohort_group
        ) %in%
          unique(
            train_metadata$cohort_group
          )
      )) {
        stop(
          "At least one CV validation cohort is absent from its training fold. ",
          "Reduce CV_FOLDS or revise the sample allocation."
        )
      }
    }
  }

  output
}

if (is.null(fold_manifest_path)) {
  fold_map <- make_stratified_fold_map(
    metadata = development_metadata,
    folds = CV_FOLDS,
    repeats = CV_REPEATS,
    seed = SEED
  )
} else {
  fold_map <- data.table::fread(fold_manifest_path)
  required_fold_fields <- c("sample_id", "Repeat", "Fold")
  if (!all(required_fold_fields %in% names(fold_map))) {
    stop("Fold manifest requires sample_id, Repeat and Fold.")
  }
  fold_map[, sample_id := as.character(sample_id)]
  fold_map[, `:=`(Repeat = as.integer(Repeat), Fold = as.integer(Fold))]
  if (anyNA(fold_map$Repeat) || anyNA(fold_map$Fold) ||
      !all(fold_map$Repeat %in% seq_len(CV_REPEATS)) ||
      !all(fold_map$Fold %in% seq_len(CV_FOLDS))) {
    stop("Fold manifest contains missing/out-of-range fold or repeat indices.")
  }
  for (r in seq_len(CV_REPEATS)) {
    current <- fold_map[Repeat == r]
    if (nrow(current) != nrow(development_metadata) ||
        anyDuplicated(current$sample_id) ||
        !setequal(current$sample_id, development_metadata$sample_id) ||
        !setequal(current$Fold, seq_len(CV_FOLDS))) {
      stop("Fold manifest must assign every development sample exactly ",
           "once in each repeat, using all CV folds.")
    }
    for (f in seq_len(CV_FOLDS)) {
      v_ids <- current[Fold == f, sample_id]
      train_md <- development_metadata[!sample_id %in% v_ids]
      valid_md <- development_metadata[sample_id %in% v_ids]
      if (nrow(valid_md) == 0L ||
          length(unique(train_md$outcome)) < 2L ||
          !all(unique(valid_md$cohort_group) %in%
               unique(train_md$cohort_group))) {
        stop("Frozen CV manifest has an invalid fold composition.")
      }
      if ("participant_key" %in% names(train_md) &&
          length(intersect(train_md$participant_key,
                           valid_md$participant_key)) > 0L) {
        stop("Frozen CV manifest splits participants between folds.")
      }
    }
  }
  # Check optional provenance columns if present.
  if ("cohort_group" %in% names(fold_map)) {
    observed_cohorts <- development_metadata$cohort_group[
      match(fold_map$sample_id, development_metadata$sample_id)]
    if (!identical(as.character(fold_map$cohort_group),
                   as.character(observed_cohorts))) {
      stop("Fold manifest cohort_group does not match metadata.")
    }
  }
  if ("outcome" %in% names(fold_map)) {
    observed_outcomes <- as.character(development_metadata$outcome[
      match(fold_map$sample_id, development_metadata$sample_id)])
    if (!identical(as.character(fold_map$outcome), observed_outcomes)) {
      stop("Fold manifest outcome does not match metadata.")
    }
  }
  # Keep output schema consistent with seeded-fold generation.
  fold_map[, cohort_group := as.character(development_metadata$cohort_group[
    match(sample_id, development_metadata$sample_id)])]
  fold_map[, outcome := as.character(development_metadata$outcome[
    match(sample_id, development_metadata$sample_id)])]
}

data.table::fwrite(
  fold_map,
  file.path(
    benchmark_dir,
    paste0(DATASET_ID, "_development_CV_fold_map.csv")
  )
)

if (VALIDATE_ONLY) {
  message("Input, split and fold checks passed for dataset ", DATASET_ID,
          ". No models fitted (--validate-only true).")
  capture.output(sessionInfo(), file.path(audit_dir,
                 paste0(DATASET_ID, "_sessionInfo_validation.txt")))
  quit(save = "no", status = 0L)
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Install these analysis packages before model fitting:\n",
       paste(missing_packages, collapse = ", "))
}

# Capture the effective dependency versions at run time, independently
# of any separately distributed renv.lock for environment restoration.
package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages,
                   function(pkg) as.character(utils::packageVersion(pkg)),
                   character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(package_versions,
  file.path(audit_dir, paste0(DATASET_ID, "_package_versions.csv")), row.names = FALSE)

# Text audit has no source-file data, but is kept private alongside other
# outputs because the run directory itself includes identifying sample data.
writeLines(c(
  paste("analysis_mode:", analysis_mode),
  paste("dataset_id:", DATASET_ID),
  paste("sampling_mode:", SAMPLING_MODE),
  paste("cases_selected:", EQUAL_CASES),
  paste("controls_selected:", EQUAL_CONTROLS),
  paste("holdout_cases:", FINAL_ICE_CASES),
  paste("holdout_controls:", FINAL_ICE_CONTROLS),
  paste("seed:", SEED),
  paste("rng_kind:", paste(RNGkind(), collapse = "/")),
  paste("cv_folds:", CV_FOLDS),
  paste("cv_repeats:", CV_REPEATS),
  paste("top_genes:", TOP_N),
  paste("cohorts:", paste(sort(unique(development_metadata$cohort_group)),
                          collapse = ",")),
  paste("manifest_supplied:", !is.null(partition_manifest_path)),
  paste("fold_map_supplied:", !is.null(fold_manifest_path)),
  paste("run_started_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(audit_dir, paste0(DATASET_ID, "_run_configuration.txt")))



# ======================================================================
# 11. VST HELPERS
# ======================================================================

calculate_poscounts_geo_means <- function(
  count_matrix
) {

  positive_number <- rowSums(
    count_matrix >
      0
  )

  if (any(
    positive_number ==
      0L
  )) {
    stop(
      "At least one retained gene is zero in every training observation."
    )
  }

  log_counts <- log(
    count_matrix
  )

  log_counts[
    !is.finite(log_counts)
  ] <- 0

  geo_means <- exp(
    rowSums(log_counts) /
      positive_number
  )

  names(geo_means) <- rownames(
    count_matrix
  )

  if (
    anyNA(geo_means) ||
      !all(is.finite(geo_means)) ||
      any(geo_means <= 0)
  ) {
    stop(
      "Invalid training geometric means."
    )
  }

  geo_means
}


fit_apply_frozen_vst <- function(
  train_counts,
  new_counts = NULL,
  train_metadata
) {

  train_counts <- train_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]

  training_total_counts <- rowSums(
    train_counts
  )

  keep_gene <- (
    is.finite(training_total_counts) &
      training_total_counts >=
        MIN_TOTAL_TRAINING_COUNT
  )

  if (sum(keep_gene) < TOP_N) {
    stop(
      "Fewer than TOP_N genes survive the training-only count filter."
    )
  }

  train_counts_filtered <- train_counts[
    keep_gene,
    ,
    drop = FALSE
  ]

  if (!is.null(new_counts)) {

    new_counts_filtered <- new_counts[
      keep_gene,
      ,
      drop = FALSE
    ]

  } else {

    new_counts_filtered <- NULL
  }

  train_coldata <- data.frame(
    outcome = factor(
      as.character(
        train_metadata$outcome
      ),
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    row.names = train_metadata$sample_id,
    check.names = FALSE
  )

  training_geo_means <- calculate_poscounts_geo_means(
    train_counts_filtered
  )

  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData = train_counts_filtered,
    colData = train_coldata,
    design = ~ outcome
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

    new_ids <- colnames(
      new_counts_filtered
    )

    new_coldata <- data.frame(
      intercept = rep(
        1L,
        length(new_ids)
      ),
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

    DESeq2::dispersionFunction(
      dds_new
    ) <- training_dispersion_function

    vst_new_object <- DESeq2::varianceStabilizingTransformation(
      dds_new,
      blind = FALSE
    )

    new_vst <- as.matrix(
      SummarizedExperiment::assay(
        vst_new_object
      )
    )

    stopifnot(
      identical(
        rownames(train_vst),
        rownames(new_vst)
      )
    )
  }

  list(
    train = train_vst,
    new = new_vst,
    state = list(
      gene_keep = keep_gene,
      genes = rownames(
        train_counts_filtered
      ),
      geo_means = training_geo_means,
      dispersion_function =
        training_dispersion_function
    )
  )
}


# ======================================================================
# 12. FROZEN VST + fSVA
# ======================================================================

fit_apply_vst_fsva <- function(
  train_counts,
  new_counts = NULL,
  train_metadata,
  seed
) {

  vst <- fit_apply_frozen_vst(
    train_counts = train_counts,
    new_counts = new_counts,
    train_metadata = train_metadata
  )

  train_vst <- vst$train
  new_vst <- vst$new

  training_variance <- matrixStats::rowVars(
    train_vst
  )

  sva_gene_keep <- (
    is.finite(training_variance) &
      training_variance >
        sqrt(
          .Machine$double.eps
        )
  )

  train_sva_input <- train_vst[
    sva_gene_keep,
    ,
    drop = FALSE
  ]

  if (!is.null(new_vst)) {
    new_sva_input <- new_vst[
      sva_gene_keep,
      ,
      drop = FALSE
    ]
  } else {
    new_sva_input <- train_sva_input
  }

  if (nrow(train_sva_input) < TOP_N) {
    stop(
      "Fewer than TOP_N variable genes available for SVA."
    )
  }

  train_coldata <- data.frame(
    outcome = factor(
      as.character(
        train_metadata$outcome
      ),
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
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
    error = function(e) {
      NA_integer_
    }
  )

  n_sv_leek <- tryCatch(
    as.integer(
      sva::num.sv(
        dat = train_sva_input,
        mod = mod_full,
        method = "leek"
      )
    ),
    error = function(e) {
      NA_integer_
    }
  )

  if (
    !is.na(n_sv_be) &&
      is.finite(n_sv_be)
  ) {

    n_sv_requested <- n_sv_be

  } else if (
    !is.na(n_sv_leek) &&
      is.finite(n_sv_leek)
  ) {

    n_sv_requested <- n_sv_leek

  } else {

    n_sv_requested <- 0L
  }

  max_sv_by_model <- max(
    0L,
    ncol(train_sva_input) -
      qr(mod_full)$rank -
      1L
  )

  max_sv_by_sample_cap <- floor(
    ncol(train_sva_input) /
      3L
  )

  max_allowed_sv <- max(
    0L,
    min(
      MAX_SV,
      max_sv_by_model,
      max_sv_by_sample_cap,
      nrow(train_sva_input) -
        1L
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

  if (n_sv_used > 0L) {

    set.seed(
      seed_for(
        seed,
        100
      )
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

    train_corrected_subset <- as.matrix(
      fsva_result$db
    )

    dimnames(
      train_corrected_subset
    ) <- dimnames(
      train_sva_input
    )

    train_fsva[
      sva_gene_keep,
      
    ] <- train_corrected_subset

    if (!is.null(new_vst)) {

      new_corrected_subset <- as.matrix(
        fsva_result$new
      )

      dimnames(
        new_corrected_subset
      ) <- dimnames(
        new_sva_input
      )

      new_fsva[
        sva_gene_keep,
        
      ] <- new_corrected_subset
    }
  }

  if (
    anyNA(train_fsva) ||
      !all(is.finite(train_fsva))
  ) {
    stop(
      "Invalid VST + fSVA training matrix."
    )
  }

  if (
    !is.null(new_fsva) &&
      (
        anyNA(new_fsva) ||
          !all(is.finite(new_fsva))
      )
  ) {
    stop(
      "Invalid VST + fSVA new-data matrix."
    )
  }

  list(
    train = train_fsva,
    new = new_fsva,
    state = c(
      vst$state,
      list(
        sva_gene_keep =
          rownames(train_vst)[
            sva_gene_keep
          ],
        nSV_BE = n_sv_be,
        nSV_Leek = n_sv_leek,
        nSV_used = n_sv_used
      )
    )
  )
}


# ======================================================================
# 13. FROZEN VST + COMBAT-REF
# ======================================================================

fit_apply_vst_combatref <- function(
  train_counts,
  new_counts = NULL,
  train_metadata
) {

  vst <- fit_apply_frozen_vst(
    train_counts = train_counts,
    new_counts = new_counts,
    train_metadata = train_metadata
  )

  expr_train <- vst$train
  expr_new <- vst$new

  train_batches <- factor(
    as.character(
      train_metadata$cohort_group
    )
  )

  cohort_counts <- sort(
    table(
      train_batches
    ),
    decreasing = TRUE
  )

  if (length(cohort_counts) < 2L) {
    stop(
      "ComBat-Ref requires at least two training cohorts."
    )
  }

  if (any(cohort_counts < 2L)) {
    stop(
      "ComBat-Ref requires at least two training observations in every cohort."
    )
  }

  ref_level <- names(
    cohort_counts
  )[[1L]]

  within_batch_residual_ss <- numeric(
    nrow(expr_train)
  )

  for (
    batch_level in unique(
      as.character(
        train_batches
      )
    )
  ) {

    batch_columns <- which(
      train_batches ==
        batch_level
    )

    batch_matrix <- expr_train[
      ,
      batch_columns,
      drop = FALSE
    ]

    batch_centred <- sweep(
      batch_matrix,
      MARGIN = 1L,
      STATS = rowMeans(
        batch_matrix
      ),
      FUN = "-"
    )

    within_batch_residual_ss <-
      within_batch_residual_ss +
      rowSums(
        batch_centred^2
      )
  }

  combat_gene_keep <- (
    is.finite(
      within_batch_residual_ss
    ) &
      within_batch_residual_ss >
        sqrt(
          .Machine$double.eps
        )
  )

  if (!any(combat_gene_keep)) {
    stop(
      "No genes contained sufficient within-cohort training variation for ComBat-Ref."
    )
  }

  expr_train_input <- expr_train[
    combat_gene_keep,
    ,
    drop = FALSE
  ]

  fit_train <- neuroCombat::neuroCombat(
    dat = expr_train_input,
    batch = as.character(
      train_batches
    ),
    mod = NULL,
    eb = TRUE,
    parametric = TRUE,
    mean.only = FALSE,
    ref.batch = ref_level,
    verbose = FALSE
  )

  adjusted_train_subset <- as.matrix(
    fit_train$dat.combat
  )

  dimnames(
    adjusted_train_subset
  ) <- dimnames(
    expr_train_input
  )

  adjusted_train <- expr_train

  adjusted_train[
    combat_gene_keep,
    
  ] <- adjusted_train_subset

  adjusted_new <- expr_new

  if (!is.null(expr_new)) {

    if (is.null(fit_train$estimates)) {
      stop(
        "neuroCombat did not return reusable training estimates."
      )
    }

    new_batches <- as.character(
      names(
        setNames(
          rep(
            NA_character_,
            ncol(expr_new)
          ),
          colnames(expr_new)
        )
      )
    )

    if (!is.null(new_counts)) {

      new_ids <- colnames(
        new_counts
      )

      if (!all(
        new_ids %in%
          colnames(expr_new)
      )) {
        stop(
          "ComBat-Ref new-data identifiers do not align."
        )
      }
    }

    stop(
      "Internal ComBat-Ref call requires new_metadata; use fit_apply_representation()."
    )
  }

  list(
    train = adjusted_train,
    new = NULL,
    state = c(
      vst$state,
      list(
        reference_cohort = ref_level,
        combat_gene_keep =
          rownames(expr_train)[
            combat_gene_keep
          ],
        combat_estimates =
          fit_train$estimates
      )
    )
  )
}


fit_apply_vst_combatref_with_new <- function(
  train_counts,
  new_counts,
  train_metadata,
  new_metadata
) {

  vst <- fit_apply_frozen_vst(
    train_counts = train_counts,
    new_counts = new_counts,
    train_metadata = train_metadata
  )

  expr_train <- vst$train
  expr_new <- vst$new

  train_batches <- as.character(
    train_metadata$cohort_group
  )

  new_batches <- as.character(
    new_metadata$cohort_group
  )

  unseen_batches <- setdiff(
    unique(new_batches),
    unique(train_batches)
  )

  if (length(unseen_batches) > 0L) {
    stop(
      "Frozen ComBat-Ref cannot transform unseen cohort(s): ",
      paste(
        unseen_batches,
        collapse = ", "
      )
    )
  }

  cohort_counts <- sort(
    table(
      train_batches
    ),
    decreasing = TRUE
  )

  if (
    length(cohort_counts) < 2L ||
      any(cohort_counts < 2L)
  ) {
    stop(
      "ComBat-Ref requires >=2 training cohorts and >=2 observations per cohort."
    )
  }

  ref_level <- names(
    cohort_counts
  )[[1L]]

  within_batch_residual_ss <- numeric(
    nrow(expr_train)
  )

  for (batch_level in unique(train_batches)) {

    batch_columns <- which(
      train_batches ==
        batch_level
    )

    batch_matrix <- expr_train[
      ,
      batch_columns,
      drop = FALSE
    ]

    batch_centred <- sweep(
      batch_matrix,
      MARGIN = 1L,
      STATS = rowMeans(
        batch_matrix
      ),
      FUN = "-"
    )

    within_batch_residual_ss <-
      within_batch_residual_ss +
      rowSums(
        batch_centred^2
      )
  }

  combat_gene_keep <- (
    is.finite(
      within_batch_residual_ss
    ) &
      within_batch_residual_ss >
        sqrt(
          .Machine$double.eps
        )
  )

  expr_train_input <- expr_train[
    combat_gene_keep,
    ,
    drop = FALSE
  ]

  expr_new_input <- expr_new[
    combat_gene_keep,
    ,
    drop = FALSE
  ]

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

  fit_new <- neuroCombat::neuroCombatFromTraining(
    dat = expr_new_input,
    batch = new_batches,
    estimates = fit_train$estimates,
    mod = NULL,
    verbose = FALSE
  )

  adjusted_train <- expr_train
  adjusted_new <- expr_new

  adjusted_train[
    combat_gene_keep,
    
  ] <- as.matrix(
    fit_train$dat.combat
  )

  adjusted_new[
    combat_gene_keep,
    
  ] <- as.matrix(
    fit_new$dat.combat
  )

  dimnames(adjusted_train) <- dimnames(
    expr_train
  )

  dimnames(adjusted_new) <- dimnames(
    expr_new
  )

  list(
    train = adjusted_train,
    new = adjusted_new,
    state = c(
      vst$state,
      list(
        reference_cohort = ref_level,
        combat_gene_keep =
          rownames(expr_train)[
            combat_gene_keep
          ],
        combat_estimates =
          fit_train$estimates
      )
    )
  )
}


# ======================================================================
# 14. FROZEN COMBAT-SEQ + LOG1P
#
# Training-estimated implementation for batches observed in training.
# No outcome is supplied to ComBat-Seq.
# ======================================================================

vec_to_matrix <- function(
  vec,
  n_columns
) {

  matrix(
    rep(
      vec,
      n_columns
    ),
    ncol = n_columns,
    byrow = FALSE
  )
}


combatseq_match_quantiles <- function(
  counts_sub,
  old_mu,
  old_phi,
  new_mu,
  new_phi
) {

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
        size = 1 /
          old_phi[
            gene_index
          ]
      )

      if (
        !is.finite(probability) ||
          abs(
            probability -
              1
          ) <
            1e-4
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
            size = 1 /
              new_phi[
                gene_index
              ]
          )
      }
    }
  }

  adjusted
}


fit_frozen_combatseq <- function(
  train_counts,
  batch_train
) {

  training_total_counts <- rowSums(
    train_counts
  )

  keep_initial <- (
    is.finite(training_total_counts) &
      training_total_counts >=
        MIN_TOTAL_TRAINING_COUNT
  )

  train_counts <- train_counts[
    keep_initial,
    ,
    drop = FALSE
  ]

  batch_train <- droplevels(
    factor(
      as.character(
        batch_train
      )
    )
  )

  if (any(
    table(batch_train) <=
      1L
  )) {
    stop(
      "Frozen ComBat-Seq requires at least two training observations in every cohort."
    )
  }

  keep_by_batch <- lapply(
    levels(batch_train),
    function(batch_level) {

      which(
        rowSums(
          train_counts[
            ,
            batch_train ==
              batch_level,
            drop = FALSE
          ]
        ) >
          0L
      )
    }
  )

  gene_keep_index <- Reduce(
    intersect,
    keep_by_batch
  )

  if (length(gene_keep_index) < TOP_N) {
    stop(
      "Fewer than TOP_N genes had non-zero counts in every training cohort for ComBat-Seq."
    )
  }

  gene_remove_index <- setdiff(
    seq_len(
      nrow(train_counts)
    ),
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

  batch_levels <- levels(
    batch_train
  )

  n_batches <- length(
    batch_levels
  )

  batch_indices <- lapply(
    batch_levels,
    function(batch_level) {
      which(
        batch_train ==
          batch_level
      )
    }
  )

  batch_sizes <- vapply(
    batch_indices,
    length,
    integer(1)
  )

  n_samples <- sum(
    batch_sizes
  )

  batch_design <- stats::model.matrix(
    ~ -1 + batch_train
  )

  colnames(
    batch_design
  ) <- batch_levels

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
        subset = nrow(
          counts_fit
        )
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
        dispersion =
          common_dispersion[
            batch_index
          ],
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
      batch_sizes[
        batch_index
      ]
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
        batch_sizes /
          n_samples,
        ncol = 1L
      )
  )

  names(alpha_g) <- rownames(
    counts_fit
  )

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
  rownames(gamma_star) <- rownames(
    counts_fit
  )

  phi_hat <- do.call(
    cbind,
    gene_dispersion
  )

  colnames(phi_hat) <- batch_levels
  rownames(phi_hat) <- rownames(
    counts_fit
  )

  phi_star <- rowMeans(
    phi_hat
  )

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
          length(
            current_columns
          )
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
      initial_genes =
        rownames(train_counts),
      fitted_genes =
        rownames(train_counts)[
          gene_keep_index
        ],
      unchanged_genes =
        rownames(train_counts)[
          gene_remove_index
        ],
      batch_levels =
        batch_levels,
      alpha_g =
        alpha_g,
      gamma_star =
        gamma_star,
      phi_hat =
        phi_hat,
      phi_star =
        phi_star
    )
  )
}


apply_frozen_combatseq <- function(
  new_counts,
  batch_new,
  state
) {

  new_counts <- new_counts[
    state$initial_genes,
    ,
    drop = FALSE
  ]

  batch_new <- as.character(
    batch_new
  )

  unseen_batches <- setdiff(
    unique(
      batch_new
    ),
    state$batch_levels
  )

  if (length(unseen_batches) > 0L) {
    stop(
      "Frozen ComBat-Seq cannot transform cohort(s) absent from training: ",
      paste(
        unseen_batches,
        collapse = ", "
      )
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

  if (any(
    !is.finite(
      new_library_offset
    )
  )) {
    stop(
      "Invalid ComBat-Seq library-size offset in new data."
    )
  }

  adjusted_fit_counts <- matrix(
    NA_real_,
    nrow = nrow(new_fit_counts),
    ncol = ncol(new_fit_counts),
    dimnames = dimnames(new_fit_counts)
  )

  for (batch_level in unique(batch_new)) {

    current_columns <- which(
      batch_new ==
        batch_level
    )

    batch_index <- match(
      batch_level,
      state$batch_levels
    )

    base_mu <- exp(
      vec_to_matrix(
        state$alpha_g,
        length(
          current_columns
        )
      ) +
        matrix(
          rep(
            new_library_offset[
              current_columns
            ],
            each =
              length(
                state$alpha_g
              )
          ),
          nrow =
            length(
              state$alpha_g
            ),
          ncol =
            length(
              current_columns
            )
        )
    )

    old_mu <- base_mu *
      exp(
        vec_to_matrix(
          state$gamma_star[
            ,
            batch_index
          ],
          length(
            current_columns
          )
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
    seq_len(
      nrow(new_counts)
    ),
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


fit_apply_frozen_combatseq_log1p <- function(
  train_counts,
  new_counts = NULL,
  train_metadata,
  new_metadata = NULL
) {

  fit <- fit_frozen_combatseq(
    train_counts = train_counts,
    batch_train =
      train_metadata$cohort_group
  )

  train_log <- log1p(
    fit$adjusted_train
  )

  new_log <- NULL

  if (!is.null(new_counts)) {

    if (is.null(new_metadata)) {
      stop(
        "new_metadata is required for frozen ComBat-Seq projection."
      )
    }

    adjusted_new <- apply_frozen_combatseq(
      new_counts = new_counts,
      batch_new =
        new_metadata$cohort_group,
      state =
        fit$state
    )

    new_log <- log1p(
      adjusted_new
    )
  }

  if (
    anyNA(train_log) ||
      !all(is.finite(train_log))
  ) {
    stop(
      "Invalid frozen ComBat-Seq + log1p training matrix."
    )
  }

  if (
    !is.null(new_log) &&
      (
        anyNA(new_log) ||
          !all(is.finite(new_log))
      )
  ) {
    stop(
      "Invalid frozen ComBat-Seq + log1p new matrix."
    )
  }

  list(
    train = train_log,
    new = new_log,
    state = fit$state
  )
}


# ======================================================================
# 15. UNIFIED REPRESENTATION WRAPPER
# ======================================================================

fit_apply_representation <- function(
  representation_id,
  train_counts,
  new_counts = NULL,
  train_metadata,
  new_metadata = NULL,
  seed
) {

  if (
    representation_id ==
      "VST_fSVA"
  ) {

    return(
      fit_apply_vst_fsva(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        seed = seed
      )
    )
  }

  if (
    representation_id ==
      "VST_ComBatRef"
  ) {

    if (is.null(new_counts)) {

      return(
        fit_apply_vst_combatref(
          train_counts = train_counts,
          new_counts = NULL,
          train_metadata = train_metadata
        )
      )
    }

    return(
      fit_apply_vst_combatref_with_new(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        new_metadata = new_metadata
      )
    )
  }

  if (
    representation_id ==
      "ComBatSeq_log1p"
  ) {

    return(
      fit_apply_frozen_combatseq_log1p(
        train_counts = train_counts,
        new_counts = new_counts,
        train_metadata = train_metadata,
        new_metadata = new_metadata
      )
    )
  }

  stop(
    "Unknown representation: ",
    representation_id
  )
}


# ======================================================================
# 16. FULL-DEVELOPMENT REPRESENTATIONS FOR DESCRIPTIVE DIAGNOSTICS
#
# Final Ice is NOT transformed here.
# ======================================================================

development_representations <- list()

representation_state_rows <- list()

for (
  representation_index in seq_along(
    representation_ids
  )
) {

  representation_id <- representation_ids[[representation_index]]

  cat(
    "\n============================================================\n",
    "FULL DEVELOPMENT DIAGNOSTIC REPRESENTATION | ",
    representation_labels[[representation_id]],
    "\n",
    "============================================================\n",
    sep = ""
  )

  current_fit <- fit_apply_representation(
    representation_id =
      representation_id,
    train_counts =
      development_counts,
    new_counts =
      NULL,
    train_metadata =
      development_metadata,
    new_metadata =
      NULL,
    seed =
      seed_for(
        200,
        representation_index
      )
  )

  development_representations[[representation_id]] <- current_fit$train

  state_row <- data.table::data.table(
    Representation_ID =
      representation_id,
    Representation =
      representation_labels[[representation_id]],
    Training_n =
      ncol(
        current_fit$train
      ),
    Genes =
      nrow(
        current_fit$train
      )
  )

  if (
    representation_id ==
      "VST_fSVA"
  ) {
    state_row[
      ,
      `:=`(
        nSV_BE =
          current_fit$state$nSV_BE,
        nSV_Leek =
          current_fit$state$nSV_Leek,
        nSV_used =
          current_fit$state$nSV_used
      )
    ]
  }

  if (
    representation_id ==
      "VST_ComBatRef"
  ) {
    state_row[
      ,
      Reference_cohort :=
        current_fit$state$reference_cohort
    ]
  }

  representation_state_rows[[representation_index]] <- state_row

  saveRDS(
    current_fit$train,
    file.path(
      diagnostic_dir,
      paste0(
        paste0(DATASET_ID, "_development_"),
        representation_id,
        ".rds"
      )
    )
  )
}

representation_state_audit <- data.table::rbindlist(
  representation_state_rows,
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  representation_state_audit,
  file.path(
    audit_dir,
    paste0(DATASET_ID, "_development_representation_state_audit.csv")
  )
)


# ======================================================================
# 17. REPRESENTATION DIAGNOSTICS: PCA + SILHOUETTE
# ======================================================================

safe_silhouette_summary <- function(
  pc_scores,
  labels,
  k
) {

  labels <- factor(
    labels
  )

  if (nlevels(labels) < 2L) {
    return(
      c(
        mean = NA_real_,
        median = NA_real_
      )
    )
  }

  scores <- pc_scores[
    ,
    seq_len(k),
    drop = FALSE
  ]

  distances <- stats::dist(
    scores
  )

  silhouette_object <- cluster::silhouette(
    as.integer(labels),
    distances
  )

  widths <- silhouette_object[
    ,
    "sil_width"
  ]

  c(
    mean = mean(
      widths
    ),
    median = stats::median(
      widths
    )
  )
}


choose_k_star <- function(
  curve
) {

  finite <- curve[
    is.finite(
      Mean_silhouette
    )
  ]

  if (nrow(finite) == 0L) {
    return(NA_integer_)
  }

  maximum <- max(
    finite$Mean_silhouette
  )

  if (maximum > 0) {

    threshold <- SILHOUETTE_NEAR_MAX *
      maximum

    candidates <- finite[
      Mean_silhouette >=
        threshold
    ]

    return(
      min(
        candidates$k
      )
    )
  }

  finite[
    which.max(
      Mean_silhouette
    ),
    k
  ]
}


representation_diagnostic_rows <- list()
silhouette_curve_rows <- list()

diagnostic_counter <- 0L
curve_counter <- 0L

for (representation_id in representation_ids) {

  expression_matrix <- development_representations[[representation_id]]

  gene_sd <- matrixStats::rowSds(
    expression_matrix
  )

  pca_gene_keep <- (
    is.finite(gene_sd) &
      gene_sd >
        sqrt(
          .Machine$double.eps
        )
  )

  pca_input <- expression_matrix[
    pca_gene_keep,
    ,
    drop = FALSE
  ]

  pca <- stats::prcomp(
    t(
      pca_input
    ),
    center = TRUE,
    scale. = TRUE
  )

  variance_percent <- (
    pca$sdev^2 /
      sum(
        pca$sdev^2
      )
  ) *
    100

  max_k <- min(
    MAX_PCA_K,
    ncol(
      pca$x
    ),
    nrow(
      pca$x
    ) -
      1L
  )

  cohort_curve_rows <- list()
  outcome_curve_rows <- list()

  for (k in 2:max_k) {

    cohort_summary <- safe_silhouette_summary(
      pc_scores = pca$x,
      labels =
        development_metadata$cohort_group,
      k = k
    )

    outcome_summary <- safe_silhouette_summary(
      pc_scores = pca$x,
      labels =
        development_metadata$outcome,
      k = k
    )

    cohort_curve_rows[[length(cohort_curve_rows) + 1L]] <- data.table::data.table(
      k = k,
      Mean_silhouette =
        cohort_summary[["mean"]],
      Median_silhouette =
        cohort_summary[["median"]]
    )

    outcome_curve_rows[[length(outcome_curve_rows) + 1L]] <- data.table::data.table(
      k = k,
      Mean_silhouette =
        outcome_summary[["mean"]],
      Median_silhouette =
        outcome_summary[["median"]]
    )
  }

  cohort_curve <- data.table::rbindlist(
    cohort_curve_rows
  )

  outcome_curve <- data.table::rbindlist(
    outcome_curve_rows
  )

  cohort_k_star <- choose_k_star(
    cohort_curve
  )

  outcome_k_star <- choose_k_star(
    outcome_curve
  )

  cohort_selected <- cohort_curve[
    k ==
      cohort_k_star
  ]

  outcome_selected <- outcome_curve[
    k ==
      outcome_k_star
  ]

  diagnostic_counter <- diagnostic_counter + 1L

  representation_diagnostic_rows[[diagnostic_counter]] <- data.table::data.table(
    Representation =
      representation_labels[[representation_id]],
    Training_n =
      nrow(
        development_metadata
      ),
    Cohorts =
      data.table::uniqueN(
        development_metadata$cohort_group
      ),
    PCA_genes =
      nrow(
        pca_input
      ),
    Cohort_k_star =
      cohort_k_star,
    Cohort_mean_silhouette =
      cohort_selected$Mean_silhouette,
    Cohort_median_silhouette =
      cohort_selected$Median_silhouette,
    Cohort_cumulative_variance_percent =
      sum(
        variance_percent[
          seq_len(
            cohort_k_star
          )
        ]
      ),
    Outcome_k_star =
      outcome_k_star,
    Outcome_mean_silhouette =
      outcome_selected$Mean_silhouette,
    Outcome_median_silhouette =
      outcome_selected$Median_silhouette,
    Outcome_cumulative_variance_percent =
      sum(
        variance_percent[
          seq_len(
            outcome_k_star
          )
        ]
      )
  )

  cohort_curve[
    ,
    `:=`(
      Representation =
        representation_labels[[representation_id]],
      Label = "Cohort"
    )
  ]

  outcome_curve[
    ,
    `:=`(
      Representation =
        representation_labels[[representation_id]],
      Label = "Outcome"
    )
  ]

  curve_counter <- curve_counter + 1L

  silhouette_curve_rows[[curve_counter]] <- data.table::rbindlist(
    list(
      cohort_curve,
      outcome_curve
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

representation_diagnostics <- data.table::rbindlist(
  representation_diagnostic_rows
)

silhouette_curves <- data.table::rbindlist(
  silhouette_curve_rows
)

data.table::fwrite(
  representation_diagnostics,
  file.path(
    diagnostic_dir,
    paste0(DATASET_ID, "_representation_diagnostics.csv")
  )
)

data.table::fwrite(
  silhouette_curves,
  file.path(
    diagnostic_dir,
    paste0(DATASET_ID, "_representation_silhouette_curves.csv")
  )
)

cat(
  "\n============================================================\n",
  "REPRESENTATION DIAGNOSTICS\n",
  "============================================================\n",
  sep = ""
)

print(
  representation_diagnostics
)


# ======================================================================
# 18. TECHNICAL / SEQUENCING-DEPTH DEPENDENCE
# ======================================================================

safe_spearman <- function(
  x,
  y
) {

  keep <- (
    is.finite(x) &
      is.finite(y)
  )

  if (
    sum(keep) <
      3L ||
      length(
        unique(
          x[
            keep
          ]
        )
      ) <
        2L ||
      length(
        unique(
          y[
            keep
          ]
        )
      ) <
        2L
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    unname(
      stats::cor(
        x[
          keep
        ],
        y[
          keep
        ],
        method = "spearman"
      )
    )
  )
}

technical_summary_rows <- list()
technical_cohort_rows <- list()

for (representation_index in seq_along(representation_ids)) {

  representation_id <- representation_ids[[representation_index]]

  transformed <- development_representations[[representation_id]]

  raw_for_representation <- development_counts[
    rownames(
      transformed
    ),
    ,
    drop = FALSE
  ]

  raw_library_size <- colSums(
    raw_for_representation,
    na.rm = TRUE
  )

  detected_genes_raw <- colSums(
    raw_for_representation >
      0,
    na.rm = TRUE
  )

  transformed_total_signal <- colSums(
    transformed,
    na.rm = TRUE
  )

  residual_signal_detection_rho <- safe_spearman(
    transformed_total_signal,
    detected_genes_raw
  )

  sequencing_depth_rho <- safe_spearman(
    raw_library_size,
    transformed_total_signal
  )

  technical_summary_rows[[representation_index]] <- data.table::data.table(
    Representation =
      representation_labels[[representation_id]],
    Residual_signal_detection_rho =
      residual_signal_detection_rho,
    Combined_sequencing_depth_rho =
      sequencing_depth_rho
  )

  qc <- data.table::data.table(
    sample_id =
      development_metadata$sample_id,
    cohort_group =
      development_metadata$cohort_group,
    raw_library_size =
      raw_library_size,
    detected_genes_raw =
      detected_genes_raw,
    transformed_total_signal =
      transformed_total_signal
  )

  current_cohort <- qc[
    ,
    .(
      n = .N,
      Sequencing_depth_rho =
        safe_spearman(
          raw_library_size,
          transformed_total_signal
        ),
      Signal_detection_rho =
        safe_spearman(
          transformed_total_signal,
          detected_genes_raw
        )
    ),
    by = cohort_group
  ]

  current_cohort[
    ,
    Representation :=
      representation_labels[[representation_id]]
  ]

  technical_cohort_rows[[representation_index]] <- current_cohort
}

technical_summary <- data.table::rbindlist(
  technical_summary_rows
)

technical_by_cohort <- data.table::rbindlist(
  technical_cohort_rows
)

data.table::fwrite(
  technical_summary,
  file.path(
    diagnostic_dir,
    paste0(DATASET_ID, "_sequencing_depth_dependence.csv")
  )
)

data.table::fwrite(
  technical_by_cohort,
  file.path(
    diagnostic_dir,
    paste0(DATASET_ID, "_sequencing_depth_by_cohort.csv")
  )
)

cat(
  "\n============================================================\n",
  "SEQUENCING-DEPTH / TECHNICAL DEPENDENCE\n",
  "============================================================\n",
  sep = ""
)

print(
  technical_summary
)

print(
  technical_by_cohort
)


# ======================================================================
# 19. MLR3 TASK + FILTER HELPERS
# ======================================================================

make_feature_task <- function(
  expression_matrix,
  metadata,
  task_id
) {

  original_genes <- rownames(
    expression_matrix
  )

  safe_genes <- make.names(
    original_genes,
    unique = TRUE
  )

  x <- data.table::as.data.table(
    t(
      expression_matrix
    )
  )

  data.table::setnames(
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
    map = data.table::data.table(
      original_gene =
        original_genes,
      safe_gene =
        safe_genes
    )
  )
}


make_filter_object <- function(
  filter_id,
  seed
) {

  if (
    filter_id ==
      "importance"
  ) {

    importance_learner <- mlr3::lrn(
      "classif.ranger",
      predict_type = "prob",
      importance = "permutation",
      num.trees = 300L,
      num.threads = 1L,
      seed = seed
    )

    return(
      mlr3filters::flt(
        "importance",
        learner =
          importance_learner
      )
    )
  }

  mlr3filters::flt(
    filter_id
  )
}


select_top_features <- function(
  train_matrix,
  validation_matrix = NULL,
  train_metadata,
  filter_id,
  top_n,
  seed
) {

  training_variance <- matrixStats::rowVars(
    train_matrix
  )

  keep <- (
    is.finite(
      training_variance
    ) &
      training_variance >
        0
  )

  train_matrix <- train_matrix[
    keep,
    ,
    drop = FALSE
  ]

  if (!is.null(validation_matrix)) {
    validation_matrix <- validation_matrix[
      keep,
      ,
      drop = FALSE
    ]
  }

  if (nrow(train_matrix) < top_n) {
    stop(
      filter_id,
      ": fewer than TOP_N variable genes."
    )
  }

  task_info <- make_feature_task(
    expression_matrix =
      train_matrix,
    metadata =
      train_metadata,
    task_id =
      paste0(
        "FS_",
        filter_id
      )
  )

  filter_object <- make_filter_object(
    filter_id = filter_id,
    seed = seed
  )

  set.seed(seed)

  filter_object$calculate(
    task_info$task
  )

  scores <- filter_object$scores

  score_table <- data.table::data.table(
    safe_gene = names(scores),
    score = as.numeric(
      scores
    )
  )

  score_table <- score_table[
    !is.na(safe_gene) &
      nzchar(safe_gene) &
      is.finite(score)
  ]

  score_table[
    ,
    original_gene :=
      task_info$map$original_gene[
        match(
          safe_gene,
          task_info$map$safe_gene
        )
      ]
  ]

  score_table <- score_table[
    !is.na(
      original_gene
    )
  ]

  data.table::setorder(
    score_table,
    -score,
    safe_gene
  )

  if (nrow(score_table) < top_n) {
    stop(
      filter_id,
      ": fewer than TOP_N finite feature scores."
    )
  }

  selected_genes <- score_table[
    seq_len(top_n),
    original_gene
  ]

  list(
    genes = selected_genes,
    scores = score_table,
    train = train_matrix[
      selected_genes,
      ,
      drop = FALSE
    ],
    validation = if (
      !is.null(
        validation_matrix
      )
    ) {
      validation_matrix[
        selected_genes,
        ,
        drop = FALSE
      ]
    } else {
      NULL
    }
  )
}


# ======================================================================
# 20. FULL-DEVELOPMENT FEATURE SETS + JACCARD CONCORDANCE
# ======================================================================

full_feature_rows <- list()
full_feature_sets <- list()
full_feature_counter <- 0L

for (representation_index in seq_along(representation_ids)) {

  representation_id <- representation_ids[[representation_index]]

  full_feature_sets[[representation_id]] <- list()

  for (fsf_index in seq_along(fsf_ids)) {

    fsf_id <- fsf_ids[[fsf_index]]

    current_fs <- select_top_features(
      train_matrix =
        development_representations[[representation_id]],
      validation_matrix =
        NULL,
      train_metadata =
        development_metadata,
      filter_id =
        fsf_id,
      top_n =
        TOP_N,
      seed =
        seed_for(
          300,
          representation_index,
          fsf_index
        )
    )

    full_feature_sets[[representation_id]][[fsf_id]] <- current_fs$genes

    full_feature_counter <-
      full_feature_counter +
      1L

    full_feature_rows[[full_feature_counter]] <- data.table::data.table(
      Representation =
        representation_labels[[representation_id]],
      FS_method =
        fsf_labels[[fsf_id]],
      Rank =
        seq_len(
          TOP_N
        ),
      Gene =
        current_fs$genes
    )
  }
}

full_feature_table <- data.table::rbindlist(
  full_feature_rows
)

data.table::fwrite(
  full_feature_table,
  file.path(
    concordance_dir,
    paste0(DATASET_ID, "_feature_sets_full_development.csv")
  )
)

jaccard_rows <- list()
jaccard_counter <- 0L

for (representation_id in representation_ids) {

  jaccard_matrix <- matrix(
    NA_real_,
    nrow = length(
      fsf_ids
    ),
    ncol = length(
      fsf_ids
    ),
    dimnames = list(
      fsf_labels[
        fsf_ids
      ],
      fsf_labels[
        fsf_ids
      ]
    )
  )

  for (i in seq_along(fsf_ids)) {

    for (j in seq_along(fsf_ids)) {

      set_a <- full_feature_sets[[representation_id]][[fsf_ids[[i]]]]

      set_b <- full_feature_sets[[representation_id]][[fsf_ids[[j]]]]

      shared <- length(
        intersect(
          set_a,
          set_b
        )
      )

      union_n <- length(
        union(
          set_a,
          set_b
        )
      )

      ji <- if (union_n > 0L) {
        shared /
          union_n
      } else {
        NA_real_
      }

      jaccard_matrix[
        i,
        j
      ] <- ji

      if (j > i) {

        jaccard_counter <-
          jaccard_counter +
          1L

        jaccard_rows[[jaccard_counter]] <- data.table::data.table(
          Representation =
            representation_labels[[representation_id]],
          FS_method_1 =
            fsf_labels[[fsf_ids[[i]]]],
          FS_method_2 =
            fsf_labels[[fsf_ids[[j]]]],
          Top_n =
            TOP_N,
          Shared_genes =
            shared,
          Jaccard_index =
            ji
        )
      }
    }
  }

  data.table::fwrite(
    data.table::as.data.table(
      jaccard_matrix,
      keep.rownames =
        "FS_method"
    ),
    file.path(
      concordance_dir,
      paste0(
        paste0(DATASET_ID, "_Jaccard_matrix_"),
        representation_id,
        ".csv"
      )
    )
  )
}

jaccard_pairwise <- data.table::rbindlist(
  jaccard_rows
)

data.table::fwrite(
  jaccard_pairwise,
  file.path(
    concordance_dir,
    paste0(DATASET_ID, "_Jaccard_pairwise.csv")
  )
)

cat(
  "\n============================================================\n",
  "FEATURE-SET JACCARD CONCORDANCE\n",
  "============================================================\n",
  sep = ""
)

print(
  jaccard_pairwise
)


# ======================================================================
# 21. LEARNERS
# ======================================================================

make_learner <- function(
  learner_id
) {

  if (
    learner_id ==
      "Naive_Bayes"
  ) {

    return(
      mlr3::lrn(
        "classif.naive_bayes",
        predict_type = "prob",
        laplace = 0
      )
    )
  }

  if (
    learner_id ==
      "Linear_SVM"
  ) {

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

  if (
    learner_id ==
      "C5_0"
  ) {

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

  if (
    learner_id ==
      "Random_Forest"
  ) {

    return(
      mlr3::lrn(
        "classif.ranger",
        predict_type = "prob",
        num.trees = 500L,
        num.threads = 1L,
        importance = "none",
        save.memory = TRUE
      )
    )
  }

  if (
    learner_id ==
      "LightGBM"
  ) {

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
        deterministic = TRUE
      )
    )
  }

  if (
    learner_id ==
      "XGBoost"
  ) {

    return(
      mlr3::lrn(
        "classif.xgboost",
        predict_type = "prob",
        nrounds = 250L,
        min_child_weight = 1,
        colsample_bytree = 0.80,
        subsample = 0.80,
        nthread = 1L,
        eval_metric = "logloss"
      )
    )
  }

  if (
    learner_id ==
      "RBF_SVM"
  ) {

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

  stop(
    "Unknown learner: ",
    learner_id
  )
}


make_classification_task <- function(
  expression_matrix,
  metadata,
  task_id
) {

  original_genes <- rownames(
    expression_matrix
  )

  safe_genes <- make.names(
    original_genes,
    unique = TRUE
  )

  x <- data.table::as.data.table(
    t(
      expression_matrix
    )
  )

  data.table::setnames(
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

  mlr3::TaskClassif$new(
    id = task_id,
    backend = x,
    target = "outcome",
    positive = POSITIVE_CLASS
  )
}


safe_prediction_metrics <- function(
  prediction
) {

  probability_case <- as.numeric(
    prediction$prob[
      ,
      POSITIVE_CLASS
    ]
  )

  truth <- as.character(
    prediction$truth
  )

  roc_object <- pROC::roc(
    response = factor(
      truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    predictor =
      probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  scores <- prediction$score(
    mlr3::msrs(
      c(
        "classif.prauc",
        "classif.bacc",
        "classif.bbrier"
      )
    )
  )

  list(
    AUROC = as.numeric(
      pROC::auc(
        roc_object
      )
    ),
    AUPRC = as.numeric(
      scores[["classif.prauc"]]
    ),
    BACC = as.numeric(
      scores[["classif.bacc"]]
    ),
    Brier = as.numeric(
      scores[["classif.bbrier"]]
    ),
    probability_case =
      probability_case
  )
}


# ======================================================================
# 22. DEVELOPMENT BENCHMARK: TRAINING-FITTED TRANSFORMATIONS
# ======================================================================

benchmark_rows <- list()
benchmark_counter <- 0L

for (
  repeat_index in seq_len(
    CV_REPEATS
  )
) {

  for (
    fold_index in seq_len(
      CV_FOLDS
    )
  ) {

    validation_ids <- fold_map[
      Repeat ==
        repeat_index &
        Fold ==
          fold_index,
      sample_id
    ]

    train_ids <- setdiff(
      development_metadata$sample_id,
      validation_ids
    )

    train_metadata <- development_metadata[
      match(
        train_ids,
        sample_id
      )
    ]

    validation_metadata <- development_metadata[
      match(
        validation_ids,
        sample_id
      )
    ]

    train_counts <- development_counts[
      ,
      train_metadata$sample_id,
      drop = FALSE
    ]

    validation_counts <- development_counts[
      ,
      validation_metadata$sample_id,
      drop = FALSE
    ]

    cat(
      "\n============================================================\n",
      "BENCHMARK | Repeat ",
      repeat_index,
      " | Fold ",
      fold_index,
      "/",
      CV_FOLDS,
      "\n",
      "============================================================\n",
      sep = ""
    )

    for (
      representation_index in seq_along(
        representation_ids
      )
    ) {

      representation_id <- representation_ids[[representation_index]]

      cat(
        "  Representation: ",
        representation_labels[[representation_id]],
        "\n",
        sep = ""
      )

      representation_attempt <- tryCatch(
        fit_apply_representation(
          representation_id =
            representation_id,
          train_counts =
            train_counts,
          new_counts =
            validation_counts,
          train_metadata =
            train_metadata,
          new_metadata =
            validation_metadata,
          seed =
            seed_for(
              400,
              repeat_index,
              fold_index,
              representation_index
            )
        ),
        error = function(e) {
          e
        }
      )

      if (
        inherits(
          representation_attempt,
          "error"
        )
      ) {

        for (fsf_id in fsf_ids) {

          for (learner_id in learner_ids) {

            benchmark_counter <-
              benchmark_counter +
              1L

            benchmark_rows[[benchmark_counter]] <- data.table::data.table(
              Repeat =
                repeat_index,
              Fold =
                fold_index,
              Representation =
                representation_labels[[representation_id]],
              Representation_ID =
                representation_id,
              FS_method =
                fsf_labels[[fsf_id]],
              FSF_ID =
                fsf_id,
              Learner =
                learner_labels[[learner_id]],
              Learner_ID =
                learner_id,
              Train_n =
                nrow(
                  train_metadata
                ),
              Validation_n =
                nrow(
                  validation_metadata
                ),
              AUROC =
                NA_real_,
              AUPRC =
                NA_real_,
              BACC =
                NA_real_,
              Brier =
                NA_real_,
              Status =
                "representation_failed",
              Error =
                conditionMessage(
                  representation_attempt
                )
            )
          }
        }

        next
      }

      for (
        fsf_index in seq_along(
          fsf_ids
        )
      ) {

        fsf_id <- fsf_ids[[fsf_index]]

        fs_attempt <- tryCatch(
          select_top_features(
            train_matrix =
              representation_attempt$train,
            validation_matrix =
              representation_attempt$new,
            train_metadata =
              train_metadata,
            filter_id =
              fsf_id,
            top_n =
              TOP_N,
            seed =
              seed_for(
                500,
                repeat_index,
                fold_index,
                representation_index,
                fsf_index
              )
          ),
          error = function(e) {
            e
          }
        )

        if (
          inherits(
            fs_attempt,
            "error"
          )
        ) {

          for (learner_id in learner_ids) {

            benchmark_counter <-
              benchmark_counter +
              1L

            benchmark_rows[[benchmark_counter]] <- data.table::data.table(
              Repeat =
                repeat_index,
              Fold =
                fold_index,
              Representation =
                representation_labels[[representation_id]],
              Representation_ID =
                representation_id,
              FS_method =
                fsf_labels[[fsf_id]],
              FSF_ID =
                fsf_id,
              Learner =
                learner_labels[[learner_id]],
              Learner_ID =
                learner_id,
              Train_n =
                nrow(
                  train_metadata
                ),
              Validation_n =
                nrow(
                  validation_metadata
                ),
              AUROC =
                NA_real_,
              AUPRC =
                NA_real_,
              BACC =
                NA_real_,
              Brier =
                NA_real_,
              Status =
                "feature_selection_failed",
              Error =
                conditionMessage(
                  fs_attempt
                )
            )
          }

          next
        }

        training_task <- make_classification_task(
          expression_matrix =
            fs_attempt$train,
          metadata =
            train_metadata,
          task_id =
            paste0(
              "train_",
              repeat_index,
              "_",
              fold_index,
              "_",
              representation_id,
              "_",
              fsf_id
            )
        )

        validation_task <- make_classification_task(
          expression_matrix =
            fs_attempt$validation,
          metadata =
            validation_metadata,
          task_id =
            paste0(
              "validation_",
              repeat_index,
              "_",
              fold_index,
              "_",
              representation_id,
              "_",
              fsf_id
            )
        )

        for (
          learner_index in seq_along(
            learner_ids
          )
        ) {

          learner_id <- learner_ids[[learner_index]]

          current_result <- tryCatch(
            {

              set.seed(
                seed_for(
                  600,
                  repeat_index,
                  fold_index,
                  representation_index,
                  fsf_index,
                  learner_index
                )
              )

              learner <- make_learner(
                learner_id
              )

              learner$train(
                training_task
              )

              prediction <- learner$predict(
                validation_task
              )

              metrics <- safe_prediction_metrics(
                prediction
              )

              data.table::data.table(
                Repeat =
                  repeat_index,
                Fold =
                  fold_index,
                Representation =
                  representation_labels[[representation_id]],
                Representation_ID =
                  representation_id,
                FS_method =
                  fsf_labels[[fsf_id]],
                FSF_ID =
                  fsf_id,
                Learner =
                  learner_labels[[learner_id]],
                Learner_ID =
                  learner_id,
                Train_n =
                  training_task$nrow,
                Validation_n =
                  validation_task$nrow,
                AUROC =
                  metrics$AUROC,
                AUPRC =
                  metrics$AUPRC,
                BACC =
                  metrics$BACC,
                Brier =
                  metrics$Brier,
                Status =
                  "success",
                Error =
                  NA_character_
              )
            },
            error = function(e) {

              data.table::data.table(
                Repeat =
                  repeat_index,
                Fold =
                  fold_index,
                Representation =
                  representation_labels[[representation_id]],
                Representation_ID =
                  representation_id,
                FS_method =
                  fsf_labels[[fsf_id]],
                FSF_ID =
                  fsf_id,
                Learner =
                  learner_labels[[learner_id]],
                Learner_ID =
                  learner_id,
                Train_n =
                  training_task$nrow,
                Validation_n =
                  validation_task$nrow,
                AUROC =
                  NA_real_,
                AUPRC =
                  NA_real_,
                BACC =
                  NA_real_,
                Brier =
                  NA_real_,
                Status =
                  "model_failed",
                Error =
                  conditionMessage(
                    e
                  )
              )
            }
          )

          benchmark_counter <-
            benchmark_counter +
            1L

          benchmark_rows[[benchmark_counter]] <- current_result
        }
      }
    }

    partial_results <- data.table::rbindlist(
      benchmark_rows,
      use.names = TRUE,
      fill = TRUE
    )

    data.table::fwrite(
      partial_results,
      file.path(
        benchmark_dir,
        paste0(DATASET_ID, "_benchmark_fold_results_PARTIAL.csv")
      )
    )
  }
}

benchmark_results <- data.table::rbindlist(
  benchmark_rows,
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  benchmark_results,
  file.path(
    benchmark_dir,
    paste0(DATASET_ID, "_benchmark_fold_results.csv")
  )
)

expected_folds <- CV_FOLDS *
  CV_REPEATS

benchmark_summary <- benchmark_results[
  Status ==
    "success" &
    is.finite(
      AUROC
    ),
  .(
    Completed_folds = .N,
    Mean_AUROC =
      mean(
        AUROC
      ),
    SD_AUROC =
      if (.N > 1L) {
        stats::sd(
          AUROC
        )
      } else {
        NA_real_
      },
    Median_AUROC =
      stats::median(
        AUROC
      ),
    Worst_AUROC =
      min(
        AUROC
      ),
    Best_AUROC =
      max(
        AUROC
      ),
    Mean_AUPRC =
      mean(
        AUPRC,
        na.rm = TRUE
      ),
    Mean_BACC =
      mean(
        BACC,
        na.rm = TRUE
      ),
    Mean_Brier =
      mean(
        Brier,
        na.rm = TRUE
      )
  ),
  by = .(
    Representation,
    Representation_ID,
    FS_method,
    FSF_ID,
    Learner,
    Learner_ID
  )
]

benchmark_summary[
  ,
  Complete :=
    Completed_folds ==
      expected_folds
]

benchmark_summary[
  Complete ==
    TRUE,
  Rank :=
    data.table::frank(
      Mean_Brier,
      ties.method = "first"
    )
]

complete_summary <- benchmark_summary[
  Complete ==
    TRUE
]

data.table::setorder(
  complete_summary,
  -Mean_AUROC,
  SD_AUROC,
  Mean_Brier
)

complete_summary[
  ,
  Rank := seq_len(.N)
]

benchmark_summary <- merge(
  benchmark_summary[
    ,
    !"Rank"
  ],
  complete_summary[
    ,
    .(
      Representation_ID,
      FSF_ID,
      Learner_ID,
      Rank
    )
  ],
  by = c(
    "Representation_ID",
    "FSF_ID",
    "Learner_ID"
  ),
  all.x = TRUE,
  sort = FALSE
)

data.table::fwrite(
  benchmark_summary,
  file.path(
    benchmark_dir,
    paste0(DATASET_ID, "_benchmark_summary.csv")
  )
)

if (nrow(complete_summary) == 0L) {
  stop(
    "No benchmark pipeline completed all required CV folds."
  )
}

locked_winner <- complete_summary[
  1L
]

data.table::fwrite(
  locked_winner,
  file.path(
    benchmark_dir,
    paste0(DATASET_ID, "_locked_winner.csv")
  )
)

cat(
  "\n============================================================\n",
  "BENCHMARK TOP 10\n",
  "============================================================\n",
  sep = ""
)

print(
  complete_summary[
    seq_len(
      min(
        10L,
        .N
      )
    )
  ]
)

cat(
  "\nLOCKED WINNER:\n"
)

print(
  locked_winner
)


# ======================================================================
# 23. CCV SPLIT TABLE
# ======================================================================

ccv_split_rows <- list()

for (target_index in seq_along(cohort_names)) {

  target_cohort <- cohort_names[[target_index]]

  source_cohorts <- setdiff(
    cohort_names,
    target_cohort
  )

  source_metadata <- development_metadata[
    cohort_group %in%
      source_cohorts
  ]

  target_metadata <- development_metadata[
    cohort_group ==
      target_cohort
  ]

  ccv_split_rows[[target_index]] <- data.table::data.table(
    Target_cohort =
      target_cohort,
    Source_cohorts =
      paste(
        source_cohorts,
        collapse = " + "
      ),
    Source_n =
      nrow(
        source_metadata
      ),
    Source_controls =
      sum(
        source_metadata$outcome ==
          NEGATIVE_CLASS
      ),
    Source_cases =
      sum(
        source_metadata$outcome ==
          POSITIVE_CLASS
      ),
    Target_n =
      nrow(
        target_metadata
      ),
    Target_controls =
      sum(
        target_metadata$outcome ==
          NEGATIVE_CLASS
      ),
    Target_cases =
      sum(
        target_metadata$outcome ==
          POSITIVE_CLASS
      )
  )
}

ccv_split_table <- data.table::rbindlist(
  ccv_split_rows
)

data.table::fwrite(
  ccv_split_table,
  file.path(
    ccv_dir,
    paste0(DATASET_ID, "_CCV_split_table.csv")
  )
)

ccv_feasibility <- data.table::data.table(
  Representation =
    representation_labels[
      representation_ids
    ],
  Strict_LOCO_estimable =
    c(
      TRUE,
      FALSE,
      FALSE
    ),
  Reason =
    c(
      paste(
        "fSVA projects nuisance structure learned from source cohorts",
        "without requiring the target cohort label during preprocessing."
      ),
      paste(
        "neuroCombatFromTraining requires the target batch level to have",
        "been represented in the training data."
      ),
      paste(
        "Frozen ComBat-Seq projection requires the target batch level to",
        "have been represented in training."
      )
    )
)

data.table::fwrite(
  ccv_feasibility,
  file.path(
    ccv_dir,
    paste0(DATASET_ID, "_CCV_representation_feasibility.csv")
  )
)

cat(
  "\n============================================================\n",
  "STRICT CCV SPLITS\n",
  "============================================================\n",
  sep = ""
)

print(
  ccv_split_table
)

cat(
  "\nCCV representation feasibility:\n"
)

print(
  ccv_feasibility
)


# ======================================================================
# 24. STRICT LOCO CCV — VST + fSVA ONLY
# ======================================================================

ccv_result_rows <- list()
ccv_result_counter <- 0L

for (target_index in seq_along(cohort_names)) {

  target_cohort <- cohort_names[[target_index]]

  source_cohorts <- setdiff(
    cohort_names,
    target_cohort
  )

  source_metadata <- development_metadata[
    cohort_group %in%
      source_cohorts
  ]

  target_metadata <- development_metadata[
    cohort_group ==
      target_cohort
  ]

  source_counts <- development_counts[
    ,
    source_metadata$sample_id,
    drop = FALSE
  ]

  target_counts <- development_counts[
    ,
    target_metadata$sample_id,
    drop = FALSE
  ]

  cat(
    "\n============================================================\n",
    "STRICT CCV | ",
    paste(
      source_cohorts,
      collapse = " + "
    ),
    " -> ",
    target_cohort,
    "\n",
    "============================================================\n",
    sep = ""
  )

  ccv_representation <- tryCatch(
    fit_apply_vst_fsva(
      train_counts =
        source_counts,
      new_counts =
        target_counts,
      train_metadata =
        source_metadata,
      seed =
        seed_for(
          700,
          target_index
        )
    ),
    error = function(e) {
      e
    }
  )

  if (
    inherits(
      ccv_representation,
      "error"
    )
  ) {

    for (fsf_id in fsf_ids) {

      for (learner_id in learner_ids) {

        ccv_result_counter <-
          ccv_result_counter +
          1L

        ccv_result_rows[[ccv_result_counter]] <- data.table::data.table(
          Target_cohort =
            target_cohort,
          Source_cohorts =
            paste(
              source_cohorts,
              collapse = " + "
            ),
          Representation =
            "VST + fSVA",
          FS_method =
            fsf_labels[[fsf_id]],
          FSF_ID =
            fsf_id,
          Learner =
            learner_labels[[learner_id]],
          Learner_ID =
            learner_id,
          Source_n =
            nrow(
              source_metadata
            ),
          Target_n =
            nrow(
              target_metadata
            ),
          AUROC =
            NA_real_,
          AUPRC =
            NA_real_,
          BACC =
            NA_real_,
          Brier =
            NA_real_,
          Status =
            "representation_failed",
          Error =
            conditionMessage(
              ccv_representation
            )
        )
      }
    }

    next
  }

  for (fsf_index in seq_along(fsf_ids)) {

    fsf_id <- fsf_ids[[fsf_index]]

    fs_attempt <- tryCatch(
      select_top_features(
        train_matrix =
          ccv_representation$train,
        validation_matrix =
          ccv_representation$new,
        train_metadata =
          source_metadata,
        filter_id =
          fsf_id,
        top_n =
          TOP_N,
        seed =
          seed_for(
            710,
            target_index,
            fsf_index
          )
      ),
      error = function(e) {
        e
      }
    )

    if (
      inherits(
        fs_attempt,
        "error"
      )
    ) {

      for (learner_id in learner_ids) {

        ccv_result_counter <-
          ccv_result_counter +
          1L

        ccv_result_rows[[ccv_result_counter]] <- data.table::data.table(
          Target_cohort =
            target_cohort,
          Source_cohorts =
            paste(
              source_cohorts,
              collapse = " + "
            ),
          Representation =
            "VST + fSVA",
          FS_method =
            fsf_labels[[fsf_id]],
          FSF_ID =
            fsf_id,
          Learner =
            learner_labels[[learner_id]],
          Learner_ID =
            learner_id,
          Source_n =
            nrow(
              source_metadata
            ),
          Target_n =
            nrow(
              target_metadata
            ),
          AUROC =
            NA_real_,
          AUPRC =
            NA_real_,
          BACC =
            NA_real_,
          Brier =
            NA_real_,
          Status =
            "feature_selection_failed",
          Error =
            conditionMessage(
              fs_attempt
            )
        )
      }

      next
    }

    source_task <- make_classification_task(
      expression_matrix =
        fs_attempt$train,
      metadata =
        source_metadata,
      task_id =
        paste0(
          "CCV_source_",
          target_index,
          "_",
          fsf_id
        )
    )

    target_task <- make_classification_task(
      expression_matrix =
        fs_attempt$validation,
      metadata =
        target_metadata,
      task_id =
        paste0(
          "CCV_target_",
          target_index,
          "_",
          fsf_id
        )
    )

    for (
      learner_index in seq_along(
        learner_ids
      )
    ) {

      learner_id <- learner_ids[[learner_index]]

      current_result <- tryCatch(
        {

          set.seed(
            seed_for(
              720,
              target_index,
              fsf_index,
              learner_index
            )
          )

          learner <- make_learner(
            learner_id
          )

          learner$train(
            source_task
          )

          prediction <- learner$predict(
            target_task
          )

          metrics <- safe_prediction_metrics(
            prediction
          )

          data.table::data.table(
            Target_cohort =
              target_cohort,
            Source_cohorts =
              paste(
                source_cohorts,
                collapse = " + "
              ),
            Representation =
              "VST + fSVA",
            FS_method =
              fsf_labels[[fsf_id]],
            FSF_ID =
              fsf_id,
            Learner =
              learner_labels[[learner_id]],
            Learner_ID =
              learner_id,
            Source_n =
              source_task$nrow,
            Target_n =
              target_task$nrow,
            AUROC =
              metrics$AUROC,
            AUPRC =
              metrics$AUPRC,
            BACC =
              metrics$BACC,
            Brier =
              metrics$Brier,
            Status =
              "success",
            Error =
              NA_character_
          )
        },
        error = function(e) {

          data.table::data.table(
            Target_cohort =
              target_cohort,
            Source_cohorts =
              paste(
                source_cohorts,
                collapse = " + "
              ),
            Representation =
              "VST + fSVA",
            FS_method =
              fsf_labels[[fsf_id]],
            FSF_ID =
              fsf_id,
            Learner =
              learner_labels[[learner_id]],
            Learner_ID =
              learner_id,
            Source_n =
              source_task$nrow,
            Target_n =
              target_task$nrow,
            AUROC =
              NA_real_,
            AUPRC =
              NA_real_,
            BACC =
              NA_real_,
            Brier =
              NA_real_,
            Status =
              "model_failed",
            Error =
              conditionMessage(
                e
              )
          )
        }
      )

      ccv_result_counter <-
        ccv_result_counter +
        1L

      ccv_result_rows[[ccv_result_counter]] <- current_result
    }
  }
}

ccv_results <- data.table::rbindlist(
  ccv_result_rows,
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  ccv_results,
  file.path(
    ccv_dir,
    paste0(DATASET_ID, "_CCV_results_VST_fSVA.csv")
  )
)

ccv_summary <- ccv_results[
  Status ==
    "success" &
    is.finite(
      AUROC
    ),
  .(
    Target_cohorts =
      .N,
    Mean_LOCO_AUROC =
      mean(
        AUROC
      ),
    SD_LOCO_AUROC =
      if (.N > 1L) {
        stats::sd(
          AUROC
        )
      } else {
        NA_real_
      },
    Median_LOCO_AUROC =
      stats::median(
        AUROC
      ),
    Worst_target_AUROC =
      min(
        AUROC
      ),
    Best_target_AUROC =
      max(
        AUROC
      ),
    Mean_LOCO_AUPRC =
      mean(
        AUPRC,
        na.rm = TRUE
      ),
    Mean_LOCO_BACC =
      mean(
        BACC,
        na.rm = TRUE
      ),
    Mean_LOCO_Brier =
      mean(
        Brier,
        na.rm = TRUE
      )
  ),
  by = .(
    FS_method,
    FSF_ID,
    Learner,
    Learner_ID
  )
]

data.table::setorder(
  ccv_summary,
  -Mean_LOCO_AUROC
)

ccv_summary[
  ,
  LOCO_rank :=
    seq_len(.N)
]

data.table::fwrite(
  ccv_summary,
  file.path(
    ccv_dir,
    paste0(DATASET_ID, "_CCV_summary_VST_fSVA.csv")
  )
)

cat(
  "\n============================================================\n",
  "STRICT CCV TOP RESULTS | VST + fSVA\n",
  "============================================================\n",
  sep = ""
)

print(
  ccv_summary[
    seq_len(
      min(
        10L,
        .N
      )
    )
  ]
)



# ======================================================================
# 25. FINAL MODEL SELECTION FROM STRICT CCV
#
# The common primary-model selection rule is:
#   - Rank pipelines using strict VST + fSVA leave-one-cohort-out CCV.
#   - Pre-specify the top four CCV pipelines BEFORE Final Ice.
#   - Candidate 1 is the primary final model.
#   - Candidates 2-4 are held-out sensitivity challengers.
#   - Final Ice is NOT used to redefine the primary model.
# ======================================================================

if (nrow(ccv_summary) < 4L) {
  stop(
    "Fewer than four successful strict-CCV pipelines are available."
  )
}

ccv_ranked_for_final <- data.table::copy(
  ccv_summary
)

data.table::setorder(
  ccv_ranked_for_final,
  -Mean_LOCO_AUROC,
  -Worst_target_AUROC,
  SD_LOCO_AUROC,
  Mean_LOCO_Brier
)

final_candidates <- ccv_ranked_for_final[
  seq_len(4L)
]

final_candidates[
  ,
  Candidate_rank := seq_len(.N)
]

final_candidates[
  ,
  `:=`(
    Representation = "VST + fSVA",
    Representation_ID = "VST_fSVA",
    Role = ifelse(
      Candidate_rank == 1L,
      "Primary CCV-selected model",
      "CCV challenger"
    )
  )
]

data.table::setcolorder(
  final_candidates,
  c(
    "Candidate_rank",
    "Role",
    "Representation",
    "Representation_ID",
    "FS_method",
    "FSF_ID",
    "Learner",
    "Learner_ID",
    setdiff(
      names(final_candidates),
      c(
        "Candidate_rank",
        "Role",
        "Representation",
        "Representation_ID",
        "FS_method",
        "FSF_ID",
        "Learner",
        "Learner_ID"
      )
    )
  )
)

data.table::fwrite(
  final_candidates,
  file.path(
    final_ice_dir,
    paste0(DATASET_ID, "_FINAL_ICE_predeclared_top4_CCV_candidates.csv")
  )
)

cat(
  "\n============================================================\n",
  "PREDECLARED FINAL-ICE CANDIDATES FROM STRICT CCV\n",
  "============================================================\n",
  sep = ""
)

print(
  final_candidates
)

cat(
  "\nPRIMARY FINAL MODEL BEFORE OPENING FINAL ICE:\n",
  final_candidates$Learner[[1L]],
  " | VST + fSVA | ",
  final_candidates$FS_method[[1L]],
  "\n",
  sep = ""
)


# ======================================================================
# 25.1 FINAL-ICE LEAKAGE AUDIT
# ======================================================================

if (
  length(
    intersect(
      development_metadata$sample_id,
      final_ice_metadata$sample_id
    )
  ) != 0L
) {
  stop(
    "LEAKAGE AUDIT FAILED: development and Final-Ice sample IDs overlap."
  )
}

if (nrow(development_metadata) != EQUAL_CASES + EQUAL_CONTROLS - FINAL_ICE_CASES - FINAL_ICE_CONTROLS) {
  stop(
    "Development sample size is inconsistent with the specified design."
  )
}

if (nrow(final_ice_metadata) != FINAL_ICE_CASES + FINAL_ICE_CONTROLS) {
  stop(
    "Holdout sample size is inconsistent with the specified design."
  )
}

if (
  sum(
    final_ice_metadata$outcome ==
      POSITIVE_CLASS
  ) != FINAL_ICE_CASES ||
    sum(
      final_ice_metadata$outcome ==
        NEGATIVE_CLASS
    ) != FINAL_ICE_CONTROLS
) {
  stop(
    "Holdout class counts do not match the configured design."
  )
}

if (
  !identical(
    colnames(development_counts),
    development_metadata$sample_id
  ) ||
    !identical(
      colnames(final_ice_counts),
      final_ice_metadata$sample_id
    )
) {
  stop(
    "Counts and metadata are not aligned before Final-Ice testing."
  )
}

message(
  "PASS: Final Ice remains separated from development."
)


# ======================================================================
# 25.2 FIT FINAL VST + fSVA USING DEVELOPMENT ONLY
# ======================================================================

set.seed(
  seed_for(
    800,
    1
  )
)

final_fsva_representation <- fit_apply_vst_fsva(
  train_counts =
    development_counts,
  new_counts =
    final_ice_counts,
  train_metadata =
    development_metadata,
  seed =
    seed_for(
      800,
      1
    )
)


# ======================================================================
# 25.3 FIT THE FOUR PREDECLARED CCV CANDIDATES
# ======================================================================

candidate_result_rows <- list()
candidate_prediction_rows <- list()
candidate_gene_rows <- list()
candidate_roc_objects <- list()

for (
  candidate_index in seq_len(
    nrow(final_candidates)
  )
) {

  current_candidate <- final_candidates[
    candidate_index
  ]

  current_fsf_id <-
    current_candidate$FSF_ID[[1L]]

  current_learner_id <-
    current_candidate$Learner_ID[[1L]]

  current_pipeline <- paste(
    current_candidate$Learner[[1L]],
    current_candidate$FS_method[[1L]],
    sep = " | "
  )

  cat(
    "\n============================================================\n",
    "FINAL ICE CANDIDATE ",
    candidate_index,
    "/4 | ",
    current_pipeline,
    "\n",
    "============================================================\n",
    sep = ""
  )

  current_features <- select_top_features(
    train_matrix =
      final_fsva_representation$train,
    validation_matrix =
      final_fsva_representation$new,
    train_metadata =
      development_metadata,
    filter_id =
      current_fsf_id,
    top_n =
      TOP_N,
    seed =
      seed_for(
        810,
        candidate_index
      )
  )

  if (length(current_features$genes) != TOP_N) {
    stop(
      "Final candidate did not return exactly TOP_N genes: ",
      current_pipeline
    )
  }

  current_train_task <- make_classification_task(
    expression_matrix =
      current_features$train,
    metadata =
      development_metadata,
    task_id =
      paste0(
        paste0(DATASET_ID, "_FINAL_TRAIN_"),
        candidate_index
      )
  )

  current_test_task <- make_classification_task(
    expression_matrix =
      current_features$validation,
    metadata =
      final_ice_metadata,
    task_id =
      paste0(
        paste0(DATASET_ID, "_FINAL_ICE_"),
        candidate_index
      )
  )

  set.seed(
    seed_for(
      820,
      candidate_index
    )
  )

  current_learner <- make_learner(
    current_learner_id
  )

  current_learner$train(
    current_train_task
  )

  current_prediction <- current_learner$predict(
    current_test_task
  )

  current_metrics <- safe_prediction_metrics(
    current_prediction
  )

  current_truth <- factor(
    as.character(
      current_prediction$truth
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  current_roc <- pROC::roc(
    response =
      current_truth,
    predictor =
      current_metrics$probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  current_auc_ci <- pROC::ci.auc(
    current_roc,
    conf.level = 0.95,
    method = "delong"
  )

  candidate_roc_objects[[current_pipeline]] <- current_roc

  candidate_result_rows[[candidate_index]] <- data.table::data.table(
    Candidate_rank =
      candidate_index,
    Role =
      current_candidate$Role[[1L]],
    Representation =
      "VST + fSVA",
    Representation_ID =
      "VST_fSVA",
    FS_method =
      current_candidate$FS_method[[1L]],
    FSF_ID =
      current_fsf_id,
    Learner =
      current_candidate$Learner[[1L]],
    Learner_ID =
      current_learner_id,
    Pipeline =
      current_pipeline,
    Train_n =
      nrow(development_metadata),
    Final_ice_n =
      nrow(final_ice_metadata),
    Final_ice_controls =
      sum(
        final_ice_metadata$outcome ==
          NEGATIVE_CLASS
      ),
    Final_ice_cases =
      sum(
        final_ice_metadata$outcome ==
          POSITIVE_CLASS
      ),
    CCV_mean_AUROC =
      current_candidate$Mean_LOCO_AUROC[[1L]],
    CCV_SD_AUROC =
      current_candidate$SD_LOCO_AUROC[[1L]],
    CCV_worst_AUROC =
      current_candidate$Worst_target_AUROC[[1L]],
    AUROC =
      current_metrics$AUROC,
    AUROC_lower95 =
      as.numeric(
        current_auc_ci[[1L]]
      ),
    AUROC_upper95 =
      as.numeric(
        current_auc_ci[[3L]]
      ),
    AUPRC =
      current_metrics$AUPRC,
    BACC =
      current_metrics$BACC,
    Brier =
      current_metrics$Brier
  )

  candidate_prediction_rows[[candidate_index]] <- data.table::data.table(
    Candidate_rank =
      candidate_index,
    Pipeline =
      current_pipeline,
    sample_id =
      final_ice_metadata$sample_id,
    cohort_group =
      final_ice_metadata$cohort_group,
    truth =
      as.character(
        current_prediction$truth
      ),
    probability_case =
      current_metrics$probability_case,
    response =
      as.character(
        current_prediction$response
      )
  )

  candidate_gene_rows[[candidate_index]] <- data.table::data.table(
    Candidate_rank =
      candidate_index,
    Pipeline =
      current_pipeline,
    Gene_rank =
      seq_len(
        length(current_features$genes)
      ),
    Gene =
      current_features$genes
  )
}

final_candidate_results <- data.table::rbindlist(
  candidate_result_rows
)

final_candidate_predictions <- data.table::rbindlist(
  candidate_prediction_rows
)

final_candidate_genes <- data.table::rbindlist(
  candidate_gene_rows
)

data.table::setorder(
  final_candidate_results,
  Candidate_rank
)

data.table::fwrite(
  final_candidate_results,
  file.path(
    final_ice_dir,
    paste0(DATASET_ID, "_FINAL_ICE_top4_CCV_candidates_results.csv")
  )
)

data.table::fwrite(
  final_candidate_predictions,
  file.path(
    final_ice_dir,
    paste0(DATASET_ID, "_FINAL_ICE_top4_CCV_candidates_predictions.csv")
  )
)

data.table::fwrite(
  final_candidate_genes,
  file.path(
    final_ice_dir,
    paste0(DATASET_ID, "_FINAL_ICE_top4_CCV_candidates_genes.csv")
  )
)

cat(
  "\n============================================================\n",
  "FINAL-ICE RESULTS | TOP FOUR STRICT-CCV CANDIDATES\n",
  "============================================================\n",
  sep = ""
)

print(
  final_candidate_results
)

cat(
  "\nPRIMARY MODEL REMAINS CANDIDATE 1, SELECTED BY STRICT CCV.\n",
  "Final-Ice AUROC is NOT used to choose a different primary model.\n",
  sep = ""
)


# ======================================================================
# 26. FIGURE OUTPUT DIRECTORY + COMMON COLOURS
# ======================================================================

figure_dir <- file.path(
  output_root,
  "06_MAIN_FIGURES"
)

permutation_dir <- file.path(
  output_root,
  "07_PERMUTATION_TESTING"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  permutation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

set1 <- RColorBrewer::brewer.pal(
  9,
  "Set1"
)

representation_order <- c(
  "VST + fSVA",
  "VST + ComBat-Ref",
  "ComBat-Seq + log1p"
)

filter_order <- c(
  "Information gain",
  "RF importance",
  "Boruta",
  "mRMR",
  "CMIM"
)

learner_order <- c(
  "Naive Bayes",
  "Linear SVM",
  "C5.0",
  "Random Forest",
  "LightGBM",
  "XGBoost",
  "RBF SVM"
)

cohort_order <- sort(
  unique(
    as.character(
      development_metadata$cohort_group
    )
  )
)

representation_colours <- stats::setNames(
  set1[1:3],
  representation_order
)

learner_colours <- stats::setNames(
  set1[
    seq_len(
      length(learner_order)
    )
  ],
  learner_order
)

cohort_colours <- stats::setNames(
  set1[
    seq_len(
      length(cohort_order)
    )
  ],
  cohort_order
)


# ======================================================================
# 26.1 FIGURE 1 — REPRESENTATION GEOMETRY
# ======================================================================

geometry_dt <- data.table::rbindlist(
  list(
    representation_diagnostics[
      ,
      .(
        Representation,
        Structure = "Cohort",
        k_star = Cohort_k_star,
        Mean_silhouette = Cohort_mean_silhouette
      )
    ],
    representation_diagnostics[
      ,
      .(
        Representation,
        Structure = "Outcome",
        k_star = Outcome_k_star,
        Mean_silhouette = Outcome_mean_silhouette
      )
    ]
  )
)

geometry_dt[
  ,
  Representation := factor(
    Representation,
    levels = representation_order
  )
]

geometry_dt[
  ,
  Structure := factor(
    Structure,
    levels = c(
      "Cohort",
      "Outcome"
    )
  )
]

geometry_dt[
  ,
  point_label := sprintf(
    "%.3f\nk*=%d",
    Mean_silhouette,
    as.integer(k_star)
  )
]

geometry_colours <- c(
  Cohort = set1[[2L]],
  Outcome = set1[[1L]]
)

p_representation_geometry <- ggplot2::ggplot(
  geometry_dt,
  ggplot2::aes(
    x = Representation,
    y = Mean_silhouette,
    group = Structure,
    colour = Structure
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55"
  ) +
  ggplot2::geom_line(
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    size = 3.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = point_label
    ),
    vjust = -0.8,
    size = 3.1,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = geometry_colours,
    name = "Structure"
  ) +
  ggplot2::labs(
    title =
      paste0("Representation geometry | ", DATASET_ID, " | ", SAMPLING_MODE),
    subtitle =
      "Outcome separation versus residual cohort clustering",
    x =
      "Transcriptomic representation",
    y =
      "Mean silhouette width"
  ) +
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold",
        colour = "black"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      ),
    panel.grid.minor =
      ggplot2::element_blank(),
    legend.position =
      "right"
  )

print(
  p_representation_geometry
)


# ======================================================================
# 26.2 FIGURE 2 — JACCARD HEATMAP
# ======================================================================

jaccard_plot_dt <- data.table::copy(
  jaccard_pairwise
)

reverse_jaccard <- jaccard_plot_dt[
  ,
  .(
    Representation,
    FS_method_1 = FS_method_2,
    FS_method_2 = FS_method_1,
    Top_n,
    Shared_genes,
    Jaccard_index
  )
]

diagonal_jaccard <- data.table::CJ(
  Representation =
    representation_order,
  FS_method_1 =
    filter_order,
  unique = TRUE
)

diagonal_jaccard[
  ,
  `:=`(
    FS_method_2 =
      FS_method_1,
    Top_n =
      TOP_N,
    Shared_genes =
      TOP_N,
    Jaccard_index =
      1
  )
]

jaccard_plot_dt <- data.table::rbindlist(
  list(
    jaccard_plot_dt,
    reverse_jaccard,
    diagonal_jaccard
  ),
  use.names = TRUE,
  fill = TRUE
)

jaccard_plot_dt <- unique(
  jaccard_plot_dt,
  by = c(
    "Representation",
    "FS_method_1",
    "FS_method_2"
  )
)

jaccard_plot_dt[
  ,
  `:=`(
    Representation =
      factor(
        Representation,
        levels = representation_order
      ),
    FS_method_1 =
      factor(
        FS_method_1,
        levels = filter_order
      ),
    FS_method_2 =
      factor(
        FS_method_2,
        levels = rev(filter_order)
      ),
    tile_label =
      sprintf(
        "%.2f",
        Jaccard_index
      )
  )
]

p_jaccard_heatmap <- ggplot2::ggplot(
  jaccard_plot_dt,
  ggplot2::aes(
    x = FS_method_1,
    y = FS_method_2,
    fill = Jaccard_index
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = tile_label
    ),
    size = 2.8,
    fontface = "bold"
  ) +
  ggplot2::facet_wrap(
    ~ Representation,
    ncol = 3
  ) +
  ggplot2::scale_fill_gradientn(
    colours =
      viridisLite::viridis(256),
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.2
      ),
    oob =
      scales::squish,
    name =
      "Jaccard\nindex"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title =
      paste0("Feature-selection concordance | ", DATASET_ID, " | ", SAMPLING_MODE),
    subtitle =
      "Pairwise Jaccard similarity among top-50 genes",
    x =
      "Feature-selection method",
    y =
      "Feature-selection method"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold",
        colour = "black"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        colour = "black"
      ),
    axis.text.y =
      ggplot2::element_text(
        colour = "black"
      ),
    panel.grid =
      ggplot2::element_blank(),
    legend.position =
      "right"
  )

print(
  p_jaccard_heatmap
)


# ======================================================================
# 26.3 FIGURE 3 — DEVELOPMENT-CV BENCHMARK HEATMAP
# ======================================================================

benchmark_heat_dt <- benchmark_summary[
  Complete ==
    TRUE &
    is.finite(
      Mean_AUROC
    )
]

benchmark_heat_dt[
  ,
  `:=`(
    Representation =
      factor(
        Representation,
        levels = representation_order
      ),
    FS_method =
      factor(
        FS_method,
        levels = filter_order
      ),
    Learner =
      factor(
        Learner,
        levels = learner_order
      ),
    tile_label =
      sprintf(
        "%.3f",
        Mean_AUROC
      )
  )
]

p_benchmark_heatmap <- ggplot2::ggplot(
  benchmark_heat_dt,
  ggplot2::aes(
    x = FS_method,
    y = Representation,
    fill = Mean_AUROC
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = tile_label
    ),
    size = 2.7,
    fontface = "bold"
  ) +
  ggplot2::facet_wrap(
    ~ Learner,
    ncol = 2,
    drop = FALSE
  ) +
  ggplot2::scale_y_discrete(
    limits = rev(
      representation_order
    ),
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(
    limits = filter_order,
    drop = FALSE
  ) +
  ggplot2::scale_fill_gradientn(
    colours =
      viridisLite::viridis(256),
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.2
      ),
    oob =
      scales::squish,
    name =
      "Mean AUROC"
  ) +
  ggplot2::labs(
    title =
      "Development CV benchmark after training-fitted preprocessing",
    subtitle =
      "Mean AUROC across cohort × outcome-stratified development CV",
    x =
      "Feature-selection method",
    y =
      "Transcriptomic representation"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold",
        colour = "black"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        colour = "black"
      ),
    axis.text.y =
      ggplot2::element_text(
        colour = "black"
      ),
    panel.grid =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        face = "bold"
      ),
    legend.position =
      "right"
  )

print(
  p_benchmark_heatmap
)


# ======================================================================
# 26.4 FIGURE 4 — INTERNAL ROBUSTNESS
# ======================================================================

robustness_dt <- benchmark_summary[
  Complete ==
    TRUE &
    is.finite(
      Mean_AUROC
    ) &
    is.finite(
      Worst_AUROC
    )
]

p_pipeline_robustness <- ggplot2::ggplot(
  robustness_dt,
  ggplot2::aes(
    x = Mean_AUROC,
    y = Worst_AUROC,
    colour = Representation,
    shape = Learner,
    size = SD_AUROC
  )
) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    colour = "grey60"
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey60"
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    colour = "grey60"
  ) +
  ggplot2::geom_point(
    alpha = 0.8
  ) +
  ggplot2::scale_colour_manual(
    values =
      representation_colours,
    name =
      "Representation"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2,
      5.5
    ),
    name =
      "Fold AUROC SD"
  ) +
  ggplot2::coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title =
      "Internal performance versus worst-fold robustness",
    subtitle =
      "Upper-right pipelines combined stronger mean discrimination with stronger worst-fold performance",
    x =
      "Mean development-CV AUROC",
    y =
      "Worst-fold AUROC"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      )
  )

print(
  p_pipeline_robustness
)


# ======================================================================
# 26.5 FIGURE 5 — STRICT CCV TOP 10
# ======================================================================

ccv_top10 <- ccv_ranked_for_final[
  seq_len(
    min(
      10L,
      .N
    )
  )
]

ccv_top10[
  ,
  Pipeline := paste(
    Learner,
    FS_method,
    sep = " | "
  )
]

ccv_top10[
  ,
  Pipeline := factor(
    Pipeline,
    levels = rev(Pipeline)
  )
]

ccv_top_prediction_dt <- merge(
  ccv_results[
    Status ==
      "success"
  ],
  ccv_top10[
    ,
    .(
      FSF_ID,
      Learner_ID,
      Pipeline
    )
  ],
  by = c(
    "FSF_ID",
    "Learner_ID"
  ),
  all = FALSE
)

ccv_top_prediction_dt[
  ,
  Pipeline := factor(
    as.character(Pipeline),
    levels = levels(
      ccv_top10$Pipeline
    )
  )
]

p_ccv_transportability <- ggplot2::ggplot(
  ccv_top_prediction_dt,
  ggplot2::aes(
    x = AUROC,
    y = Pipeline,
    colour = Target_cohort
  )
) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    colour = "grey55"
  ) +
  ggplot2::geom_point(
    position =
      ggplot2::position_dodge(
        width = 0.5
      ),
    size = 3
  ) +
  ggplot2::scale_colour_manual(
    values =
      cohort_colours,
    name =
      "Held-out cohort"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.1
      )
  ) +
  ggplot2::labs(
    title =
      paste0("Strict cross-cohort transportability | ", DATASET_ID),
    subtitle =
      "Top 10 pipelines ranked by leave-one-cohort-out AUROC; VST + fSVA only",
    x =
      "Held-out cohort AUROC",
    y =
      "Feature-selection method | classifier"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      ),
    panel.grid.major.y =
      ggplot2::element_blank()
  )

print(
  p_ccv_transportability
)


# ======================================================================
# 26.6 FIGURE 6 — DEVELOPMENT CV VERSUS STRICT CCV
# ======================================================================

cv_fsva <- benchmark_summary[
  Complete ==
    TRUE &
    Representation_ID ==
      "VST_fSVA",
  .(
    FS_method,
    FSF_ID,
    Learner,
    Learner_ID,
    CV_Mean_AUROC =
      Mean_AUROC,
    CV_SD_AUROC =
      SD_AUROC,
    CV_Worst_AUROC =
      Worst_AUROC
  )
]

ccv_for_merge <- ccv_summary[
  ,
  .(
    FS_method,
    FSF_ID,
    Learner,
    Learner_ID,
    CCV_Mean_AUROC =
      Mean_LOCO_AUROC,
    CCV_SD_AUROC =
      SD_LOCO_AUROC,
    CCV_Worst_AUROC =
      Worst_target_AUROC
  )
]

cv_ccv_dt <- merge(
  cv_fsva,
  ccv_for_merge,
  by = c(
    "FS_method",
    "FSF_ID",
    "Learner",
    "Learner_ID"
  ),
  all = FALSE
)

cv_ccv_dt[
  ,
  Delta_CCV_minus_CV :=
    CCV_Mean_AUROC -
    CV_Mean_AUROC
]

p_cv_vs_ccv <- ggplot2::ggplot(
  cv_ccv_dt,
  ggplot2::aes(
    x = CV_Mean_AUROC,
    y = CCV_Mean_AUROC,
    colour = Learner,
    shape = FS_method
  )
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45"
  ) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    colour = "grey65"
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    colour = "grey65"
  ) +
  ggplot2::geom_point(
    size = 3.5,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(
    values =
      learner_colours,
    name =
      "Classifier"
  ) +
  ggplot2::coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title =
      "Internal learnability versus strict cohort transportability",
    subtitle =
      "VST + fSVA pipelines only; points below the diagonal lost performance under unseen-cohort transfer",
    x =
      "Mean development-CV AUROC",
    y =
      "Mean strict CCV AUROC",
    shape =
      "Feature selection"
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      )
  )

print(
  p_cv_vs_ccv
)


# ======================================================================
# 27. FIGURE 7 — FINAL-ICE ROC CURVES + 95% BOOTSTRAP ROC BANDS
#
# AUROC confidence interval:
#   DeLong 95% CI
#
# Shaded ROC band:
#   stratified bootstrap 95% CI
# ======================================================================

ROC_BOOTSTRAPS <- 500L

specificity_grid <- seq(
  0,
  1,
  length.out = 101L
)

roc_curve_rows <- list()
roc_ribbon_rows <- list()

for (
  candidate_index in seq_len(
    nrow(final_candidates)
  )
) {

  current_result <- final_candidate_results[
    Candidate_rank ==
      candidate_index
  ]

  current_pipeline <-
    current_result$Pipeline[[1L]]

  current_roc <-
    candidate_roc_objects[[current_pipeline]]

  current_legend <- paste0(
    "Rank ",
    candidate_index,
    ": ",
    current_pipeline,
    "\nAUROC ",
    sprintf(
      "%.3f",
      current_result$AUROC[[1L]]
    ),
    " (95% DeLong CI ",
    sprintf(
      "%.3f",
      current_result$AUROC_lower95[[1L]]
    ),
    "–",
    sprintf(
      "%.3f",
      current_result$AUROC_upper95[[1L]]
    ),
    ")"
  )

  roc_curve_rows[[candidate_index]] <- data.table::data.table(
    Candidate_rank =
      candidate_index,
    Pipeline =
      current_pipeline,
    Legend =
      current_legend,
    FPR =
      1 -
      current_roc$specificities,
    Sensitivity =
      current_roc$sensitivities
  )

  set.seed(
    seed_for(
      880,
      candidate_index
    )
  )

  current_ci_se <- pROC::ci.se(
    current_roc,
    specificities =
      specificity_grid,
    conf.level =
      0.95,
    boot.n =
      ROC_BOOTSTRAPS,
    boot.stratified =
      TRUE
  )

  current_ci_matrix <- as.matrix(
    current_ci_se
  )

  roc_ribbon_rows[[candidate_index]] <- data.table::data.table(
    Candidate_rank =
      candidate_index,
    Pipeline =
      current_pipeline,
    Legend =
      current_legend,
    Specificity =
      as.numeric(
        rownames(
          current_ci_matrix
        )
      ),
    Sensitivity_lower =
      current_ci_matrix[
        ,
        1L
      ],
    Sensitivity_median =
      current_ci_matrix[
        ,
        2L
      ],
    Sensitivity_upper =
      current_ci_matrix[
        ,
        3L
      ]
  )

  roc_ribbon_rows[[candidate_index]][
    ,
    FPR :=
      1 -
      Specificity
  ]
}

roc_curve_dt <- data.table::rbindlist(
  roc_curve_rows
)

roc_ribbon_dt <- data.table::rbindlist(
  roc_ribbon_rows
)

pipeline_names <- final_candidate_results[
  order(Candidate_rank),
  Pipeline
]

pipeline_colours <- stats::setNames(
  set1[
    1:4
  ],
  pipeline_names
)

legend_labels <- unique(
  roc_curve_dt[
    order(Candidate_rank),
    .(
      Pipeline,
      Legend
    )
  ]
)

legend_lookup <- stats::setNames(
  legend_labels$Legend,
  legend_labels$Pipeline
)

p_final_ice_roc <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data =
      roc_ribbon_dt,
    ggplot2::aes(
      x = FPR,
      ymin = Sensitivity_lower,
      ymax = Sensitivity_upper,
      fill = Pipeline,
      group = Pipeline
    ),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45"
  ) +
  ggplot2::geom_step(
    data =
      roc_curve_dt,
    ggplot2::aes(
      x = FPR,
      y = Sensitivity,
      colour = Pipeline,
      group = Pipeline
    ),
    linewidth = 1.15,
    direction = "hv"
  ) +
  ggplot2::scale_colour_manual(
    values =
      pipeline_colours,
    breaks =
      names(pipeline_colours),
    labels =
      legend_lookup[
        names(pipeline_colours)
      ],
    name = NULL
  ) +
  ggplot2::scale_fill_manual(
    values =
      pipeline_colours,
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.2
      ),
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.2
      ),
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title =
      paste0("Held-out performance of strict-CCV-selected ", DATASET_ID, " pipelines"),
    subtitle =
      sprintf("Final holdout: n = %d (%d cases, %d controls); shading shows 95%% bootstrap ROC bands", nrow(final_ice_metadata), FINAL_ICE_CASES, FINAL_ICE_CONTROLS),
    x =
      "False-positive rate (1 - specificity)",
    y =
      "True-positive rate (sensitivity)"
  ) +
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      ),
    panel.grid.minor =
      ggplot2::element_blank(),
    panel.border =
      ggplot2::element_rect(
        fill = NA,
        colour = "black"
      ),
    legend.position =
      "right"
  )

print(
  p_final_ice_roc
)


# ======================================================================
# 28. QUICK FINAL-ICE PERMUTATION TEST
#
# 1000 cohort-stratified label permutations.
# Prediction scores remain fixed.
# ======================================================================

FINAL_PERMUTATIONS <- 1000L

calculate_fixed_auroc <- function(
  truth,
  probability
) {

  truth_factor <- factor(
    truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  if (
    length(
      unique(truth_factor)
    ) < 2L
  ) {
    return(
      NA_real_
    )
  }

  current_roc <- pROC::roc(
    response =
      truth_factor,
    predictor =
      probability,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(
    pROC::auc(
      current_roc
    )
  )
}

permutation_template <- final_candidate_predictions[
  Candidate_rank ==
    1L,
  .(
    sample_id,
    cohort_group,
    truth
  )
]

data.table::setorder(
  permutation_template,
  sample_id
)

prediction_matrix <- matrix(
  NA_real_,
  nrow =
    nrow(permutation_template),
  ncol =
    nrow(final_candidates)
)

colnames(
  prediction_matrix
) <- final_candidate_results[
  order(Candidate_rank),
  Pipeline
]

for (
  candidate_index in seq_len(
    nrow(final_candidates)
  )
) {

  current_prediction <- final_candidate_predictions[
    Candidate_rank ==
      candidate_index,
    .(
      sample_id,
      probability_case
    )
  ]

  data.table::setorder(
    current_prediction,
    sample_id
  )

  if (
    !identical(
      current_prediction$sample_id,
      permutation_template$sample_id
    )
  ) {
    stop(
      "Final-Ice prediction sample alignment failed."
    )
  }

  prediction_matrix[
    ,
    candidate_index
  ] <- current_prediction$probability_case
}

observed_auroc <- vapply(
  seq_len(
    ncol(prediction_matrix)
  ),
  function(candidate_index) {

    calculate_fixed_auroc(
      truth =
        permutation_template$truth,
      probability =
        prediction_matrix[
          ,
          candidate_index
        ]
    )
  },
  FUN.VALUE =
    numeric(1L)
)

null_auroc <- matrix(
  NA_real_,
  nrow =
    FINAL_PERMUTATIONS,
  ncol =
    ncol(prediction_matrix)
)

colnames(
  null_auroc
) <- colnames(
  prediction_matrix
)

set.seed(
  seed_for(
    900,
    1
  )
)

for (
  permutation_index in seq_len(
    FINAL_PERMUTATIONS
  )
) {

  permuted_truth <-
    permutation_template$truth

  for (
    cohort_name in unique(
      permutation_template$cohort_group
    )
  ) {

    current_rows <- which(
      permutation_template$cohort_group ==
        cohort_name
    )

    permuted_truth[
      current_rows
    ] <- sample(
      permutation_template$truth[
        current_rows
      ],
      size =
        length(current_rows),
      replace = FALSE
    )
  }

  for (
    candidate_index in seq_len(
      ncol(prediction_matrix)
    )
  ) {

    null_auroc[
      permutation_index,
      candidate_index
    ] <- calculate_fixed_auroc(
      truth =
        permuted_truth,
      probability =
        prediction_matrix[
          ,
          candidate_index
        ]
    )
  }

  if (
    permutation_index %% 250L ==
      0L ||
      permutation_index ==
        FINAL_PERMUTATIONS
  ) {
    message(
      "Final-Ice permutation ",
      permutation_index,
      "/",
      FINAL_PERMUTATIONS
    )
  }
}

permutation_results <- data.table::data.table(
  Candidate_rank =
    final_candidate_results$Candidate_rank,
  Pipeline =
    final_candidate_results$Pipeline,
  Observed_AUROC =
    observed_auroc,
  Null_mean =
    colMeans(
      null_auroc,
      na.rm = TRUE
    ),
  Null_median =
    apply(
      null_auroc,
      2L,
      stats::median,
      na.rm = TRUE
    ),
  Null_lower95 =
    apply(
      null_auroc,
      2L,
      stats::quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
  Null_upper95 =
    apply(
      null_auroc,
      2L,
      stats::quantile,
      probs = 0.975,
      na.rm = TRUE
    )
)

permutation_results[
  ,
  Empirical_p := vapply(
    seq_len(.N),
    function(candidate_index) {

      (
        1 +
          sum(
            null_auroc[
              ,
              candidate_index
            ] >=
              Observed_AUROC[[candidate_index]],
            na.rm = TRUE
          )
      ) /
        (
          FINAL_PERMUTATIONS +
            1
        )
    },
    FUN.VALUE =
      numeric(1L)
  )
]

permutation_results[
  ,
  BH_adjusted_p :=
    stats::p.adjust(
      Empirical_p,
      method = "BH"
    )
]

data.table::fwrite(
  permutation_results,
  file.path(
    permutation_dir,
    paste0(DATASET_ID, "_FINAL_ICE_top4_CCV_permutation_summary.csv")
  )
)

null_plot_dt <- data.table::rbindlist(
  lapply(
    seq_len(
      ncol(null_auroc)
    ),
    function(candidate_index) {

      data.table::data.table(
        Candidate_rank =
          candidate_index,
        Pipeline =
          colnames(
            null_auroc
          )[[candidate_index]],
        Null_AUROC =
          null_auroc[
            ,
            candidate_index
          ]
      )
    }
  )
)

null_plot_dt[
  ,
  Pipeline := factor(
    Pipeline,
    levels =
      final_candidate_results[
        order(Candidate_rank),
        Pipeline
      ]
  )
]

permutation_results[
  ,
  Pipeline := factor(
    Pipeline,
    levels =
      levels(
        null_plot_dt$Pipeline
      )
  )
]

permutation_results[
  ,
  permutation_label := paste0(
    "Observed = ",
    sprintf(
      "%.3f",
      Observed_AUROC
    ),
    "\np = ",
    format.pval(
      Empirical_p,
      digits = 2,
      eps = 0.001
    ),
    "; BH = ",
    format.pval(
      BH_adjusted_p,
      digits = 2,
      eps = 0.001
    )
  )
]

p_final_permutation <- ggplot2::ggplot(
  null_plot_dt,
  ggplot2::aes(
    x = Null_AUROC
  )
) +
  ggplot2::geom_histogram(
    bins = 25,
    fill = "grey82",
    colour = "white"
  ) +
  ggplot2::geom_vline(
    data =
      permutation_results,
    ggplot2::aes(
      xintercept =
        Observed_AUROC,
      colour =
        Pipeline
    ),
    linewidth = 1.15
  ) +
  ggplot2::geom_text(
    data =
      permutation_results,
    ggplot2::aes(
      x =
        Observed_AUROC,
      y =
        Inf,
      label =
        permutation_label,
      colour =
        Pipeline
    ),
    hjust = -0.04,
    vjust = 1.2,
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(
    ~ Pipeline,
    ncol = 2
  ) +
  ggplot2::scale_colour_manual(
    values =
      pipeline_colours,
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks =
      seq(
        0,
        1,
        by = 0.1
      )
  ) +
  ggplot2::labs(
    title =
      "Held-out discrimination relative to the permutation null",
    subtitle =
      "1000 cohort-stratified permutations; prediction scores remained frozen",
    x =
      "AUROC under permuted outcome labels",
    y =
      "Permutation count"
  ) +
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    axis.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text =
      ggplot2::element_text(
        colour = "black"
      ),
    strip.text =
      ggplot2::element_text(
        face = "bold"
      ),
    panel.grid.minor =
      ggplot2::element_blank()
  )

print(
  p_final_permutation
)


# ======================================================================
# 29. PAIRED DELONG COMPARISONS AMONG FOUR FINAL-ICE CANDIDATES
# ======================================================================

candidate_pairs <- utils::combn(
  seq_len(
    nrow(final_candidates)
  ),
  2L
)

delong_rows <- list()

for (
  pair_index in seq_len(
    ncol(candidate_pairs)
  )
) {

  candidate_a <-
    candidate_pairs[
      1L,
      pair_index
    ]

  candidate_b <-
    candidate_pairs[
      2L,
      pair_index
    ]

  pipeline_a <- final_candidate_results[
    Candidate_rank ==
      candidate_a,
    Pipeline
  ]

  pipeline_b <- final_candidate_results[
    Candidate_rank ==
      candidate_b,
    Pipeline
  ]

  roc_a <-
    candidate_roc_objects[[pipeline_a]]

  roc_b <-
    candidate_roc_objects[[pipeline_b]]

  current_test <- pROC::roc.test(
    roc_a,
    roc_b,
    method = "delong",
    paired = TRUE
  )

  delong_rows[[pair_index]] <- data.table::data.table(
    Pipeline_A =
      pipeline_a,
    Pipeline_B =
      pipeline_b,
    AUROC_A =
      as.numeric(
        pROC::auc(roc_a)
      ),
    AUROC_B =
      as.numeric(
        pROC::auc(roc_b)
      ),
    Delta_A_minus_B =
      as.numeric(
        pROC::auc(roc_a)
      ) -
      as.numeric(
        pROC::auc(roc_b)
      ),
    DeLong_p =
      current_test$p.value
  )
}

delong_results <- data.table::rbindlist(
  delong_rows
)

delong_results[
  ,
  BH_adjusted_p :=
    stats::p.adjust(
      DeLong_p,
      method = "BH"
    )
]

data.table::fwrite(
  delong_results,
  file.path(
    final_ice_dir,
    paste0(DATASET_ID, "_FINAL_ICE_top4_CCV_paired_DeLong.csv")
  )
)

print(
  delong_results
)


# ======================================================================
# 30. SAVE ALL MAIN FIGURES
# ======================================================================

pdf_device <- if (
  capabilities(
    "cairo"
  )
) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

figure_manifest <- list(
  Figure_1 =
    p_representation_geometry,
  Figure_2 =
    p_jaccard_heatmap,
  Figure_3 =
    p_benchmark_heatmap,
  Figure_4 =
    p_pipeline_robustness,
  Figure_5 =
    p_ccv_transportability,
  Figure_6 =
    p_cv_vs_ccv,
  Figure_7 =
    p_final_ice_roc,
  Figure_8 =
    p_final_permutation
)

figure_sizes <- list(
  Figure_1 =
    c(
      205,
      145
    ),
  Figure_2 =
    c(
      260,
      135
    ),
  Figure_3 =
    c(
      225,
      295
    ),
  Figure_4 =
    c(
      215,
      175
    ),
  Figure_5 =
    c(
      225,
      190
    ),
  Figure_6 =
    c(
      205,
      175
    ),
  Figure_7 =
    c(
      245,
      180
    ),
  Figure_8 =
    c(
      225,
      180
    )
)

figure_file_stems <- c(
  Figure_1 =
    paste0("Figure_1_", DATASET_ID, "_", SAMPLING_MODE, "_representation_geometry"),
  Figure_2 =
    paste0("Figure_2_", DATASET_ID, "_", SAMPLING_MODE, "_Jaccard_heatmap"),
  Figure_3 =
    paste0("Figure_3_", DATASET_ID, "_", SAMPLING_MODE, "_development_CV_heatmap"),
  Figure_4 =
    paste0("Figure_4_", DATASET_ID, "_", SAMPLING_MODE, "_pipeline_robustness"),
  Figure_5 =
    paste0("Figure_5_", DATASET_ID, "_", SAMPLING_MODE, "_strict_CCV_transportability"),
  Figure_6 =
    paste0("Figure_6_", DATASET_ID, "_", SAMPLING_MODE, "_CV_vs_CCV"),
  Figure_7 =
    paste0("Figure_7_", DATASET_ID, "_", SAMPLING_MODE, "_top4_CCV_Final_Ice_ROC"),
  Figure_8 =
    paste0("Figure_8_", DATASET_ID, "_", SAMPLING_MODE, "_Final_Ice_permutation")
)

for (
  figure_name in names(
    figure_manifest
  )
) {

  current_plot <-
    figure_manifest[[figure_name]]

  current_size <-
    figure_sizes[[figure_name]]

  current_stem <-
    figure_file_stems[[figure_name]]

  ggplot2::ggsave(
    filename =
      file.path(
        figure_dir,
        paste0(
          current_stem,
          ".png"
        )
      ),
    plot =
      current_plot,
    width =
      current_size[[1L]],
    height =
      current_size[[2L]],
    units =
      "mm",
    dpi =
      600,
    bg =
      "white"
  )

  ggplot2::ggsave(
    filename =
      file.path(
        figure_dir,
        paste0(
          current_stem,
          ".tiff"
        )
      ),
    plot =
      current_plot,
    width =
      current_size[[1L]],
    height =
      current_size[[2L]],
    units =
      "mm",
    dpi =
      600,
    compression =
      "lzw",
    bg =
      "white"
  )

  ggplot2::ggsave(
    filename =
      file.path(
        figure_dir,
        paste0(
          current_stem,
          ".pdf"
        )
      ),
    plot =
      current_plot,
    width =
      current_size[[1L]],
    height =
      current_size[[2L]],
    units =
      "mm",
    device =
      pdf_device,
    bg =
      "white"
  )
}


# ======================================================================
# 31. FINAL AUDITS
# ======================================================================

failure_summary <- data.table::rbindlist(
  list(
    benchmark_results[
      Status !=
        "success",
      .(
        Stage =
          "Development benchmark",
        Status,
        Error,
        N =
          .N
      ),
      by = .(
        Status,
        Error
      )
    ],
    ccv_results[
      Status !=
        "success",
      .(
        Stage =
          "Strict CCV",
        Status,
        Error,
        N =
          .N
      ),
      by = .(
        Status,
        Error
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

data.table::fwrite(
  failure_summary,
  file.path(
    audit_dir,
    paste0(DATASET_ID, "_failure_summary.csv")
  )
)

writeLines(
  c(
    paste0("TRACE-Omics | ", DATASET_ID, " | ", SAMPLING_MODE, " sampling"),
    "",
    "DESIGN",
    sprintf("Aligned %s cohort data -> selected %d cases + %d controls.", DATASET_ID, EQUAL_CASES, EQUAL_CONTROLS),
    sprintf("Development: %d cases/%d controls; holdout: %d cases/%d controls.", EQUAL_CASES-FINAL_ICE_CASES, EQUAL_CONTROLS-FINAL_ICE_CONTROLS, FINAL_ICE_CASES, FINAL_ICE_CONTROLS),
    "",
    "REPRESENTATIONS",
    "VST + fSVA",
    "VST + ComBat-Ref",
    "ComBat-Seq + log1p",
    "",
    "FEATURE SELECTION",
    "Information gain, RF importance, Boruta, mRMR, CMIM; top 50; fold-internal.",
    "",
    "BENCHMARK",
    paste0(
      CV_FOLDS,
      "-fold CV x ",
      CV_REPEATS,
      " repeat(s); fixed untuned NB/Linear-SVM/C5.0/RF/LGBM/XGB/RBF-SVM."
    ),
    "",
    "STRICT CCV",
    "Strict leave-one-cohort-out CCV is run only for VST + fSVA.",
    "",
    "FINAL MODEL SELECTION",
    "The top four strict-CCV pipelines are predeclared before Final Ice.",
    "Strict-CCV rank 1 remains the primary model regardless of Final-Ice ranking.",
    "",
    "FINAL ICE",
    sprintf("The same %d holdout samples are evaluated by four predeclared candidates.", FINAL_ICE_CASES+FINAL_ICE_CONTROLS),
    "AUROC CIs are DeLong; ROC ribbons are stratified bootstrap confidence bands.",
    "",
    "PERMUTATION",
    "1000 cohort-stratified label permutations with frozen Final-Ice prediction scores.",
    "",
    "REPRODUCIBILITY AND DATA GOVERNANCE",
    paste("Run status:", analysis_mode),
    "Exact historical reproduction requires the original frozen partition and fold maps.",
    "Input fingerprints and selected sample IDs are saved in this RESTRICTED directory.",
    "Do not upload this output directory, its split files or predictions to GitHub/Zenodo.",
    "Only independently reviewed aggregate figures/tables may be made public.",
    "This script has not been independently re-executed or scientifically validated by the editor.",
    "",
    "IMPORTANT",
    sprintf("Selected-sample case prevalence: %.3f (mode: %s).", EQUAL_CASES/(EQUAL_CASES+EQUAL_CONTROLS), SAMPLING_MODE),
    "AUPRC, Brier and calibration reflect the sampled/held-out case prevalence, not necessarily deployment prevalence."
  ),
  con =
    file.path(
      output_root,
      paste0("README_", DATASET_ID, "_", SAMPLING_MODE, "_EXPERIMENT.txt")
    )
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      audit_dir,
      paste0(DATASET_ID, "_sessionInfo.txt")
    )
)


# ======================================================================
# 32. FINAL CONSOLE SUMMARY
# ======================================================================

cat(
  "\n\n============================================================\n",
  paste0(DATASET_ID, " ", SAMPLING_MODE, " EXPERIMENT COMPLETE\n"),
  "============================================================\n",
  sep = ""
)

cat(
  "\nTRAIN / TEST SPLIT\n\n"
)

print(
  split_summary
)

cat(
  "\nREPRESENTATION DIAGNOSTICS\n\n"
)

print(
  representation_diagnostics
)

cat(
  "\nSEQUENCING-DEPTH DIAGNOSTICS\n\n"
)

print(
  technical_summary
)

cat(
  "\nDEVELOPMENT-CV WINNER — DESCRIPTIVE ONLY\n\n"
)

print(
  locked_winner
)

cat(
  "\nSTRICT CCV TOP FOUR — PREDECLARED FINAL CANDIDATES\n\n"
)

print(
  final_candidates
)

cat(
  "\nFINAL ICE TOP-FOUR RESULTS\n\n"
)

print(
  final_candidate_results
)

cat(
  "\nFINAL ICE PERMUTATION RESULTS\n\n"
)

print(
  permutation_results
)

cat(
  "\nPRIMARY FINAL MODEL:\n",
  final_candidates$Learner[[1L]],
  " | VST + fSVA | ",
  final_candidates$FS_method[[1L]],
  "\n",
  sep = ""
)

cat(
  "\nOUTPUT DIRECTORY:\n",
  output_root,
  "\n",
  sep = ""
)

if (nrow(failure_summary) == 0L) {

  cat(
    "\nPASS: no benchmark or CCV failures were recorded.\n"
  )

} else {

  cat(
    "\nWARNING: some fits failed. Inspect:\n",
    file.path(
      audit_dir,
      paste0(DATASET_ID, "_failure_summary.csv")
    ),
    "\n",
    sep = ""
  )
}

cat(
  "\nFINAL ICE WAS USED ONLY AFTER THE TOP FOUR STRICT-CCV PIPELINES ",
  "WERE PREDECLARED.\n",
  sep = ""
)
