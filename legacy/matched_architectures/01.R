# ======================================================================
# 0. LOAD CD DATA ONLY
#
# This script NEVER reuses development_counts/development_metadata from
# another disease. Any existing generic objects are removed first, then the
# frozen CD 80/20 split is loaded explicitly from:
#
#   CD_FULL_DATA_FINAL_MODEL/01_split
# ======================================================================

DISEASE <- "CD"

for (
  object_name in c(
    "development_counts",
    "development_metadata",
    "final_ice_counts",
    "final_ice_metadata"
  )
) {

  if (
    exists(
      object_name,
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {

    rm(
      list = object_name,
      envir = .GlobalEnv
    )
  }
}

gc(
  verbose = FALSE
)


USER_HOME <- Sys.getenv(
  "USERPROFILE",
  unset = path.expand("~")
)


PROJECT_DIR_FROM_DISK <- file.path(
  USER_HOME,
  "Desktop",
  paste0(
    DISEASE,
    "_FULL_DATA_FINAL_MODEL"
  )
)


SPLIT_DIR_FROM_DISK <- file.path(
  PROJECT_DIR_FROM_DISK,
  "01_split"
)


if (
  !dir.exists(
    SPLIT_DIR_FROM_DISK
  )
) {

  stop(
    "The CD split directory does not exist:\n",
    SPLIT_DIR_FROM_DISK
  )
}


resolve_disease_file <- function(
  candidates,
  label
) {

  paths <- file.path(
    SPLIT_DIR_FROM_DISK,
    candidates
  )


  hits <- paths[
    file.exists(
      paths
    )
  ]


  if (
    length(
      hits
    ) !=
      1L
  ) {

    stop(
      "Could not uniquely resolve ",
      label,
      " for CD.\n",
      "Expected one of:\n",
      paste(
        paths,
        collapse = "\n"
      ),
      "\nFound ",
      length(
        hits
      ),
      " matching file(s)."
    )
  }


  normalizePath(
    hits[[1L]],
    winslash = "/",
    mustWork = TRUE
  )
}


development_counts_path <- resolve_disease_file(
  candidates = c(
    "CD_development_raw_counts.rds"
  ),
  label = "development raw counts"
)


final_ice_counts_path <- resolve_disease_file(
  candidates = c(
    "CD_FINAL_ICE_raw_counts.rds"
  ),
  label = "FINAL ICE raw counts"
)


development_metadata_path <- resolve_disease_file(
  candidates = c(
    "CD_development_metadata.csv",
    "CD_development_metadata.rds"
  ),
  label = "development metadata"
)


final_ice_metadata_path <- resolve_disease_file(
  candidates = c(
    "CD_FINAL_ICE_metadata.csv",
    "CD_FINAL_ICE_metadata.rds"
  ),
  label = "FINAL ICE metadata"
)


development_counts <- readRDS(
  development_counts_path
)


final_ice_counts <- readRDS(
  final_ice_counts_path
)


read_metadata_file <- function(
  path
) {

  if (
    grepl(
      "\\.rds$",
      path,
      ignore.case = TRUE
    )
  ) {

    return(
      readRDS(
        path
      )
    )
  }


  data.table::fread(
    path
  )
}


development_metadata <- read_metadata_file(
  development_metadata_path
)


final_ice_metadata <- read_metadata_file(
  final_ice_metadata_path
)


# ----------------------------------------------------------------------
# HARD SOURCE AUDIT
# ----------------------------------------------------------------------

expected_project_fragment <- paste0(
  "/",
  DISEASE,
  "_FULL_DATA_FINAL_MODEL/"
)


all_source_paths <- gsub(
  "\\\\",
  "/",
  c(
    development_counts_path,
    final_ice_counts_path,
    development_metadata_path,
    final_ice_metadata_path
  )
)


if (
  any(
    !grepl(
      expected_project_fragment,
      all_source_paths,
      fixed = TRUE
    )
  )
) {

  stop(
    "SOURCE AUDIT FAILED: at least one input is not from the CD project."
  )
}


if (
  "disease" %in%
    names(
      development_metadata
    )
) {

  observed_disease <- unique(
    toupper(
      trimws(
        as.character(
          development_metadata$disease
        )
      )
    )
  )


  observed_disease <- observed_disease[
    nzchar(
      observed_disease
    )
  ]


  if (
    length(
      observed_disease
    ) >
      0L &&
      any(
        observed_disease !=
          DISEASE
      )
  ) {

    stop(
      "Development metadata disease audit failed. Observed: ",
      paste(
        observed_disease,
        collapse = ", "
      ),
      "; expected: ",
      DISEASE,
      "."
    )
  }
}


if (
  "disease" %in%
    names(
      final_ice_metadata
    )
) {

  observed_ice_disease <- unique(
    toupper(
      trimws(
        as.character(
          final_ice_metadata$disease
        )
      )
    )
  )


  observed_ice_disease <- observed_ice_disease[
    nzchar(
      observed_ice_disease
    )
  ]


  if (
    length(
      observed_ice_disease
    ) >
      0L &&
      any(
        observed_ice_disease !=
          DISEASE
      )
  ) {

    stop(
      "FINAL ICE metadata disease audit failed. Observed: ",
      paste(
        observed_ice_disease,
        collapse = ", "
      ),
      "; expected: ",
      DISEASE,
      "."
    )
  }
}


cat(
  "\n============================================================\n",
  "CD DATA SOURCE AUDIT PASSED\n",
  "============================================================\n",
  "Development counts:\n  ",
  development_counts_path,
  "\nDevelopment metadata:\n  ",
  development_metadata_path,
  "\nFINAL ICE counts:\n  ",
  final_ice_counts_path,
  "\nFINAL ICE metadata:\n  ",
  final_ice_metadata_path,
  "\nDevelopment n: ",
  ncol(
    development_counts
  ),
  "\nFINAL ICE n: ",
  ncol(
    final_ice_counts
  ),
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# CD LIGHTGBM vs TABPFN FOUNDATION-MODEL EXPERIMENT
# ======================================================================
#
# Loads the frozen CD development/ICE split directly from disk.
#
# MODELLING OBJECTS LOADED ABOVE:
#   development_counts
#   development_metadata
#   final_ice_counts
#   final_ice_metadata
#
# DESIGN:
#   A. Shared 5-fold nested CV on development data
#      - joint cohort x outcome stratification
#      - log2(CPM + 1)
#      - training-only RF permutation feature selection
#      - same 100 genes supplied to LightGBM and TabPFN
#      - LightGBM tuned in 3-fold inner CV
#      - TabPFN used without task-specific hyperparameter tuning
#
#   B. Strict leave-one-cohort-out (LOCO) transportability experiment
#      - entire target cohort excluded from preprocessing, feature selection,
#        LightGBM tuning, and fitting
#      - same source-selected genes supplied to both models
#
#   C. Diagnostic final ICE comparison
#      - final LightGBM hyperparameters chosen from nested inner-CV only
#      - feature selection performed on development data only
#      - ICE used once for diagnostic evaluation
#
# OUTPUT:
#   CD_FULL_DATA_FINAL_MODEL/103_LIGHTGBM_VS_TABPFN
#
# IMPORTANT:
#   Final Python models are saved with reticulate::py_save_object(), not
#   saveRDS(), so they can be reloaded in a new R session.
# ======================================================================


# ======================================================================
# 0. DO NOT CLEAR THE CURRENT R SESSION
# ======================================================================

gc(verbose = FALSE)


# ======================================================================
# 1. R PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "matrixStats",
  "ranger",
  "pROC",
  "PRROC",
  "reticulate",
  "ggplot2",
  "patchwork",
  "RColorBrewer",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these R packages first:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(PRROC)
  library(reticulate)
  library(patchwork)
})


# ======================================================================
# 2. REQUIRED CURRENT-SESSION OBJECTS
# ======================================================================

required_objects <- c(
  "development_counts",
  "development_metadata",
  "final_ice_counts",
  "final_ice_metadata"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    FUN.VALUE = logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Missing required current-session objects:\n",
    paste(missing_objects, collapse = "\n")
  )
}


# ======================================================================
# 3. SETTINGS
# ======================================================================

SEED_FOUNDATION <- 7832L

POSITIVE_CLASS <- if (exists("POSITIVE_CLASS", inherits = TRUE)) {
  get("POSITIVE_CLASS", inherits = TRUE)
} else {
  "case"
}

NEGATIVE_CLASS <- if (exists("NEGATIVE_CLASS", inherits = TRUE)) {
  get("NEGATIVE_CLASS", inherits = TRUE)
} else {
  "control"
}

OUTER_FOLDS <- 5L
INNER_FOLDS <- 3L
FOUNDATION_TOP_N <- 100L
RF_PREFILTER <- 3000L
RF_IMPORTANCE_TREES <- 250L

MODEL_ROC_BOOTSTRAPS <- 2000L
PR_BOOTSTRAPS <- 2000L

ROC_SPECIFICITY_GRID <- seq(0, 1, length.out = 101L)
PR_RECALL_GRID <- seq(0, 1, length.out = 101L)

PROJECT_DIR <- file.path(
  Sys.getenv("USERPROFILE"),
  "Desktop",
  "CD_FULL_DATA_FINAL_MODEL"
)

FOUNDATION_DIR <- file.path(
  PROJECT_DIR,
  "103_LIGHTGBM_VS_TABPFN"
)

PUBLICATION_DIR <- file.path(
  FOUNDATION_DIR,
  "publication"
)

GENERALISATION_DIR <- file.path(
  FOUNDATION_DIR,
  "strict_LOCO_generalisation"
)

ICE_COMPARISON_DIR <- file.path(
  FOUNDATION_DIR,
  "diagnostic_FINAL_ICE"
)

MODEL_DIR <- file.path(
  FOUNDATION_DIR,
  "saved_models"
)

for (current_dir in c(
  FOUNDATION_DIR,
  PUBLICATION_DIR,
  GENERALISATION_DIR,
  ICE_COMPARISON_DIR,
  MODEL_DIR
)) {
  dir.create(
    current_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

SET1 <- RColorBrewer::brewer.pal(9, "Set1")

FOUNDATION_LABEL <- "TabPFN"

MODEL_ORDER <- c(
  "LightGBM",
  FOUNDATION_LABEL
)

MODEL_COLOURS <- stats::setNames(
  SET1[1:2],
  MODEL_ORDER
)


# ======================================================================
# 4. STANDARDISE AND AUDIT METADATA
# ======================================================================

prepare_foundation_metadata <- function(metadata) {

  metadata <- data.table::as.data.table(
    data.table::copy(metadata)
  )

  if (!"sample_id" %in% names(metadata)) {
    if (!"geo_accession" %in% names(metadata)) {
      stop("Metadata contains neither sample_id nor geo_accession.")
    }

    metadata[
      ,
      sample_id := as.character(geo_accession)
    ]
  }

  if (!"cohort_group" %in% names(metadata)) {
    if (!"cohort" %in% names(metadata)) {
      stop("Metadata contains neither cohort_group nor cohort.")
    }

    metadata[
      ,
      cohort_group := as.character(cohort)
    ]
  }

  metadata[
    ,
    sample_id := trimws(as.character(sample_id))
  ]

  metadata[
    ,
    cohort_group := trimws(as.character(cohort_group))
  ]

  metadata[
    ,
    outcome := factor(
      tolower(trimws(as.character(outcome))),
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    )
  ]

  if (anyNA(metadata$outcome)) {
    stop("Invalid outcome labels detected.")
  }

  if (anyDuplicated(metadata$sample_id) > 0L) {
    stop("Duplicated sample IDs detected.")
  }

  metadata
}


development_metadata <- prepare_foundation_metadata(
  development_metadata
)

final_ice_metadata <- prepare_foundation_metadata(
  final_ice_metadata
)


development_counts <- as.matrix(development_counts)
final_ice_counts <- as.matrix(final_ice_counts)

storage.mode(development_counts) <- "integer"
storage.mode(final_ice_counts) <- "integer"

if (!identical(
  rownames(development_counts),
  rownames(final_ice_counts)
)) {
  if (all(
    rownames(development_counts) %in%
      rownames(final_ice_counts)
  )) {
    final_ice_counts <- final_ice_counts[
      rownames(development_counts),
      ,
      drop = FALSE
    ]
  } else {
    stop("Development and ICE gene sets cannot be aligned.")
  }
}

if (!all(
  development_metadata$sample_id %in%
    colnames(development_counts)
)) {
  stop("Some development sample IDs are missing from development_counts.")
}

if (!all(
  final_ice_metadata$sample_id %in%
    colnames(final_ice_counts)
)) {
  stop("Some ICE sample IDs are missing from final_ice_counts.")
}

development_counts <- development_counts[
  ,
  development_metadata$sample_id,
  drop = FALSE
]

final_ice_counts <- final_ice_counts[
  ,
  final_ice_metadata$sample_id,
  drop = FALSE
]

cohort_audit <- development_metadata[
  ,
  .(
    n = .N,
    Cases = sum(outcome == POSITIVE_CLASS),
    Controls = sum(outcome == NEGATIVE_CLASS)
  ),
  by = cohort_group
]

cohort_audit[
  ,
  Both_classes := Cases > 0L & Controls > 0L
]

print(cohort_audit)

if (any(!cohort_audit$Both_classes)) {
  stop(
    "Every CD development cohort must contain both classes for strict LOCO AUROC."
  )
}


# ======================================================================
# HARD DISEASE IDENTITY CHECK BEFORE ANY MODEL FITTING
# ======================================================================

if (
  !all(
    development_metadata$sample_id %in%
      colnames(
        development_counts
      )
  ) ||
    !all(
      final_ice_metadata$sample_id %in%
        colnames(
          final_ice_counts
        )
    )
) {

  stop(
    "CD sample-ID audit failed before modelling."
  )
}


cat(
  "\nPASS: modelling objects are the frozen CD development/ICE data.\n",
  "Development n = ",
  nrow(
    development_metadata
  ),
  "; ICE n = ",
  nrow(
    final_ice_metadata
  ),
  "\n",
  sep = ""
)

COHORT_ORDER <- sort(
  unique(development_metadata$cohort_group)
)

cat(
  "\n============================================================\n",
  "CD FOUNDATION-MODEL DATA AUDIT\n",
  "============================================================\n",
  "Development samples: ", nrow(development_metadata), "\n",
  "Diagnostic ICE samples: ", nrow(final_ice_metadata), "\n",
  "Genes: ", nrow(development_counts), "\n",
  "Development cohorts: ", paste(COHORT_ORDER, collapse = ", "), "\n",
  "============================================================\n",
  sep = ""
)

fwrite(
  cohort_audit,
  file.path(
    FOUNDATION_DIR,
    "00_development_cohort_audit.csv"
  )
)


# ======================================================================
# 5. DEDICATED CD PYTHON BACKEND
#
# CD uses its OWN Conda environment:
#   cd_foundation_models_clean
#
# Run each disease experiment in a fresh RStudio session because reticulate
# cannot switch Python interpreters after Python has been initialised.
# ======================================================================

FOUNDATION_CONDA_ENV <- "cd_foundation_models_clean"
FOUNDATION_PYTHON_VERSION <- "3.10"
FOUNDATION_PYTHON_PACKAGES <- c(
  "numpy==1.22.4",
  "scikit-learn==1.2.2",
  "lightgbm==4.6.0",
  "torch==2.4.1",
  "tabpfn==0.1.11"
)

if (reticulate::py_available(initialize = FALSE)) {
  active_config <- reticulate::py_config()
  active_python <- normalizePath(active_config$python, winslash = "/", mustWork = TRUE)
  if (!grepl(paste0("/envs/", FOUNDATION_CONDA_ENV, "/"), tolower(active_python), fixed = TRUE)) {
    stop(
      "Python is already initialised to the wrong environment.\n\n",
      "Active Python:\n  ", active_python, "\n\n",
      "Required environment:\n  ", FOUNDATION_CONDA_ENV, "\n\n",
      "Close ALL RStudio windows, reopen RStudio, and source this CD script first."
    )
  }
}

conda_environment_table <- reticulate::conda_list()
if (!FOUNDATION_CONDA_ENV %in% conda_environment_table$name) {
  if (reticulate::py_available(initialize = FALSE)) {
    stop("The dedicated CD environment does not exist, but Python is already initialised. Restart RStudio and rerun.")
  }
  cat("\nCreating dedicated CD Python 3.10 environment: ", FOUNDATION_CONDA_ENV, "\n", sep = "")
  reticulate::conda_create(envname = FOUNDATION_CONDA_ENV, python_version = FOUNDATION_PYTHON_VERSION)
  cat("\nInstalling pinned packages into ", FOUNDATION_CONDA_ENV, "...\n", sep = "")
  reticulate::py_install(packages = FOUNDATION_PYTHON_PACKAGES, envname = FOUNDATION_CONDA_ENV, method = "conda", pip = TRUE)
}

FOUNDATION_PYTHON <- reticulate::conda_python(FOUNDATION_CONDA_ENV)
if (is.null(FOUNDATION_PYTHON) || !file.exists(FOUNDATION_PYTHON)) stop("Could not resolve Python for ", FOUNDATION_CONDA_ENV, ".")
requested_python <- normalizePath(FOUNDATION_PYTHON, winslash = "/", mustWork = TRUE)
if (!grepl(paste0("/envs/", FOUNDATION_CONDA_ENV, "/"), tolower(requested_python), fixed = TRUE)) stop("Resolved Python does not belong to the dedicated CD environment:\n", requested_python)

TABPFN_CHECKPOINT <- file.path(dirname(FOUNDATION_PYTHON), "Lib", "site-packages", "tabpfn", "models_diff", "prior_diff_real_checkpoint_n_0_epoch_42.cpkt")
EXPECTED_TABPFN_CHECKPOINT_SIZE <- 103350223L
checkpoint_is_valid <- file.exists(TABPFN_CHECKPOINT) && isTRUE(file.info(TABPFN_CHECKPOINT)$size == EXPECTED_TABPFN_CHECKPOINT_SIZE)
if (!checkpoint_is_valid) {
  if (file.exists(TABPFN_CHECKPOINT)) {
    bad_checkpoint <- paste0(TABPFN_CHECKPOINT, ".corrupt")
    if (file.exists(bad_checkpoint)) unlink(bad_checkpoint)
    file.rename(TABPFN_CHECKPOINT, bad_checkpoint)
  }
  dir.create(dirname(TABPFN_CHECKPOINT), recursive = TRUE, showWarnings = FALSE)
  cat("\nDownloading the verified legacy TabPFN checkpoint for CD (~103 MB)...\n")
  utils::download.file(url = "https://huggingface.co/spaces/TabPFN/TabPFNPrediction/resolve/b392eb0f0cfb5efe64af2be58b0a60521abed8f8/TabPFN/models_diff/prior_diff_real_checkpoint_n_0_epoch_42.cpkt?download=true", destfile = TABPFN_CHECKPOINT, mode = "wb", method = "libcurl", quiet = FALSE)
}
checkpoint_size_observed <- file.info(TABPFN_CHECKPOINT)$size
if (!isTRUE(checkpoint_size_observed == EXPECTED_TABPFN_CHECKPOINT_SIZE)) stop("Invalid TabPFN checkpoint in ", FOUNDATION_CONDA_ENV, ".\nObserved size: ", checkpoint_size_observed, "\nExpected size: ", EXPECTED_TABPFN_CHECKPOINT_SIZE)
message("PASS: CD TabPFN checkpoint is present with the expected size.")

if (!reticulate::py_available(initialize = FALSE)) {
  Sys.setenv(RETICULATE_PYTHON = FOUNDATION_PYTHON)
  reticulate::use_python(FOUNDATION_PYTHON, required = TRUE)
}
python_config <- reticulate::py_config()
active_python <- normalizePath(python_config$python, winslash = "/", mustWork = TRUE)
if (!identical(tolower(active_python), tolower(requested_python))) stop("reticulate selected the wrong interpreter.\nSelected:\n  ", active_python, "\nRequired:\n  ", requested_python)

python_smoke_code <- paste0(
  "import sys; import numpy; import sklearn; import lightgbm; import torch; import tabpfn; ",
  "print('Python:', sys.version); print('NumPy:', numpy.__version__); ",
  "print('scikit-learn:', sklearn.__version__); print('LightGBM:', lightgbm.__version__); ",
  "print('Torch:', torch.__version__); print('TabPFN: available')"
)
smoke_output <- system2(command = FOUNDATION_PYTHON, args = c("-c", shQuote(python_smoke_code)), stdout = TRUE, stderr = TRUE)
smoke_status <- attr(smoke_output, "status")
cat("\nCD Python external smoke test:\n", paste(smoke_output, collapse = "\n"), "\n", sep = "")
if (!is.null(smoke_status) && smoke_status != 0L) stop("The dedicated CD foundation-model environment failed its smoke test.\n", paste(smoke_output, collapse = "\n"))

np <- reticulate::import("numpy", convert = FALSE)
lightgbm_py <- reticulate::import("lightgbm", convert = FALSE)
tabpfn_py <- reticulate::import("tabpfn", convert = FALSE)
torch_py <- reticulate::import("torch", convert = TRUE)
TABPFN_DEVICE <- if (isTRUE(torch_py$cuda$is_available())) "cuda" else "cpu"
cat("\n============================================================\n", "CD FOUNDATION-MODEL PYTHON BACKEND\n", "============================================================\n", "Environment: ", FOUNDATION_CONDA_ENV, "\n", "Python: ", python_config$python, "\n", "Version: ", as.character(python_config$version), "\n", "TabPFN device: ", TABPFN_DEVICE, "\n", "============================================================\n", sep = "")
message("PASS: dedicated CD Python environment is active.")


# ======================================================================
# 6. REPRESENTATION: SAMPLE-WISE log2(CPM + 1)
# ======================================================================

foundation_logcpm_transform <- function(counts) {

  library_size <- colSums(counts)

  if (any(library_size <= 0)) {
    stop("Invalid RNA-seq library size detected.")
  }

  cpm <- sweep(
    counts,
    MARGIN = 2L,
    STATS = library_size,
    FUN = "/"
  ) * 1e6

  log2(cpm + 1)
}


foundation_fit_apply_logcpm <- function(
  train_counts,
  new_counts = NULL,
  train_metadata
) {

  train_counts <- train_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]

  if (!is.null(new_counts)) {
    new_counts <- new_counts[
      rownames(train_counts),
      ,
      drop = FALSE
    ]
  }

  training_total_counts <- rowSums(train_counts)

  keep_gene <- is.finite(training_total_counts) &
    training_total_counts >= 10L

  train_matrix <- foundation_logcpm_transform(
    train_counts
  )[
    keep_gene,
    ,
    drop = FALSE
  ]

  new_matrix <- NULL

  if (!is.null(new_counts)) {
    new_matrix <- foundation_logcpm_transform(
      new_counts
    )[
      keep_gene,
      ,
      drop = FALSE
    ]
  }

  list(
    train = train_matrix,
    new = new_matrix,
    genes = rownames(train_matrix)
  )
}


# ======================================================================
# 7. TRAINING-ONLY RF PERMUTATION FEATURE SELECTION
# ======================================================================

foundation_class_weights <- function(outcome) {

  y <- factor(
    as.character(outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  counts <- table(y)

  weights <- length(y) /
    (2 * as.numeric(counts))

  names(weights) <- names(counts)

  weights
}


select_foundation_rf_features <- function(
  train_matrix,
  train_metadata,
  top_n,
  seed
) {

  gene_variance <- matrixStats::rowVars(
    train_matrix
  )

  valid <- is.finite(gene_variance) &
    gene_variance > 0

  train_matrix <- train_matrix[
    valid,
    ,
    drop = FALSE
  ]

  gene_variance <- gene_variance[
    valid
  ]

  if (nrow(train_matrix) > RF_PREFILTER) {
    variance_order <- order(
      gene_variance,
      decreasing = TRUE
    )

    keep <- variance_order[
      seq_len(RF_PREFILTER)
    ]

    train_matrix <- train_matrix[
      keep,
      ,
      drop = FALSE
    ]
  }

  x <- as.data.frame(
    t(train_matrix),
    check.names = FALSE
  )

  y <- factor(
    as.character(train_metadata$outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  class_weights <- foundation_class_weights(y)

  set.seed(seed)

  rf <- ranger::ranger(
    x = x,
    y = y,
    probability = TRUE,
    importance = "permutation",
    num.trees = RF_IMPORTANCE_TREES,
    mtry = max(
      1L,
      floor(sqrt(ncol(x)))
    ),
    class.weights = class_weights,
    num.threads = 0L,
    seed = seed
  )

  importance_values <- rf$variable.importance

  importance_values[
    !is.finite(importance_values)
  ] <- -Inf

  selected_genes <- names(
    head(
      sort(
        importance_values,
        decreasing = TRUE
      ),
      top_n
    )
  )

  if (length(selected_genes) != top_n) {
    stop(
      "RF feature selection returned ",
      length(selected_genes),
      " genes; expected ",
      top_n,
      "."
    )
  }

  selected_genes
}


# ======================================================================
# 8. JOINT COHORT x OUTCOME FOLD GENERATOR
# ======================================================================

make_joint_stratified_folds <- function(
  metadata,
  number_of_folds,
  seed
) {

  metadata <- data.table::copy(metadata)

  metadata[
    ,
    stratum := paste(
      cohort_group,
      outcome,
      sep = "__"
    )
  ]

  fold_id <- integer(
    nrow(metadata)
  )

  set.seed(seed)

  for (current_stratum in unique(metadata$stratum)) {

    current_index <- which(
      metadata$stratum == current_stratum
    )

    current_index <- sample(
      current_index
    )

    assigned_folds <- rep(
      seq_len(number_of_folds),
      length.out = length(current_index)
    )

    fold_id[
      current_index
    ] <- assigned_folds
  }

  fold_id
}


# ======================================================================
# 9. LIGHTGBM GRID
# ======================================================================

lightgbm_grid <- data.table::data.table(

  Config_ID = paste0(
    "LGBM_",
    sprintf("%02d", 1:10)
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

print(lightgbm_grid)


# ======================================================================
# 10. R -> NUMPY
# ======================================================================

to_numpy_x <- function(x) {
  np$array(
    x,
    dtype = "float32"
  )
}


to_numpy_y <- function(y) {

  y_numeric <- as.integer(
    as.character(y) ==
      POSITIVE_CLASS
  )

  np$array(
    y_numeric,
    dtype = "int64"
  )
}


# ======================================================================
# 11. PYTHON MODEL FIT / PREDICT
# ======================================================================

fit_lightgbm_object <- function(
  x_train,
  y_train,
  parameters,
  seed
) {

  model <- lightgbm_py$LGBMClassifier(
    objective = "binary",
    n_estimators = as.integer(parameters$n_estimators[[1L]]),
    learning_rate = parameters$learning_rate[[1L]],
    num_leaves = as.integer(parameters$num_leaves[[1L]]),
    max_depth = as.integer(parameters$max_depth[[1L]]),
    min_child_samples = as.integer(parameters$min_child_samples[[1L]]),
    subsample = parameters$subsample[[1L]],
    subsample_freq = 1L,
    colsample_bytree = parameters$colsample_bytree[[1L]],
    reg_alpha = parameters$reg_alpha[[1L]],
    reg_lambda = parameters$reg_lambda[[1L]],
    class_weight = "balanced",
    random_state = as.integer(seed),
    n_jobs = -1L,
    verbosity = -1L
  )

  model$fit(
    to_numpy_x(x_train),
    to_numpy_y(y_train)
  )

  model
}


predict_python_classifier <- function(
  model,
  x_test
) {

  probability_py <- model$predict_proba(
    to_numpy_x(x_test)
  )

  probability_matrix <- reticulate::py_to_r(
    probability_py
  )

  if (!is.matrix(probability_matrix)) {
    probability_matrix <- as.matrix(
      probability_matrix
    )
  }

  if (ncol(probability_matrix) != 2L) {
    stop(
      "Expected two probability columns; received ",
      ncol(probability_matrix),
      "."
    )
  }

  probability_case <- as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )

  if (any(!is.finite(probability_case))) {
    stop("Model returned non-finite probabilities.")
  }

  probability_case
}


fit_predict_lightgbm <- function(
  x_train,
  y_train,
  x_test,
  parameters,
  seed
) {

  model <- fit_lightgbm_object(
    x_train = x_train,
    y_train = y_train,
    parameters = parameters,
    seed = seed
  )

  predict_python_classifier(
    model = model,
    x_test = x_test
  )
}


fit_tabpfn_object <- function(
  x_train,
  y_train,
  seed
) {

  seed_integer <- as.integer(
    seed %%
      .Machine$integer.max
  )

  np$random$seed(
    seed_integer
  )

  model <- tabpfn_py$TabPFNClassifier(
    device = TABPFN_DEVICE,
    seed = seed_integer
  )

  model$fit(
    to_numpy_x(x_train),
    to_numpy_y(y_train)
  )

  model
}


fit_predict_tabpfn <- function(
  x_train,
  y_train,
  x_test,
  seed
) {

  model <- fit_tabpfn_object(
    x_train = x_train,
    y_train = y_train,
    seed = seed
  )

  predict_python_classifier(
    model = model,
    x_test = x_test
  )
}


# ======================================================================
# 12. PERFORMANCE METRICS
# ======================================================================

evaluate_binary_predictions <- function(
  truth,
  probability,
  threshold = 0.5
) {

  truth <- factor(
    as.character(truth),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  probability <- as.numeric(
    probability
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

  pr_object <- PRROC::pr.curve(
    scores.class0 = probability[
      truth == POSITIVE_CLASS
    ],
    scores.class1 = probability[
      truth == NEGATIVE_CLASS
    ],
    curve = FALSE
  )

  response <- ifelse(
    probability >= threshold,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  )

  sensitivity <- mean(
    response[
      truth == POSITIVE_CLASS
    ] == POSITIVE_CLASS
  )

  specificity <- mean(
    response[
      truth == NEGATIVE_CLASS
    ] == NEGATIVE_CLASS
  )

  data.table::data.table(
    AUROC = as.numeric(
      pROC::auc(roc_object)
    ),
    AUPRC = pr_object$auc.integral,
    BACC = mean(
      c(
        sensitivity,
        specificity
      )
    ),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Brier = mean(
      (
        probability -
          as.integer(
            truth == POSITIVE_CLASS
          )
      )^2
    )
  )
}


# ======================================================================
# 13. PUBLICATION HELPERS
# ======================================================================

publication_theme <- function() {

  ggplot2::theme_bw(
    base_size = 11
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        colour = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9.5,
        colour = "black"
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        colour = "black"
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        colour = "black"
      ),
      legend.text = ggplot2::element_text(
        colour = "black"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black"
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = "black"
      )
    )
}


make_calibration_data <- function(
  truth,
  probability,
  label,
  number_of_bins = 8L
) {

  calibration_dt <- data.table::data.table(
    truth = as.integer(
      truth == POSITIVE_CLASS
    ),
    probability = as.numeric(
      probability
    )
  )

  calibration_dt[
    ,
    Probability_rank := data.table::frank(
      probability,
      ties.method = "average"
    )
  ]

  calibration_dt[
    ,
    Bin := ceiling(
      Probability_rank /
        .N *
        number_of_bins
    )
  ]

  calibration_dt[
    Bin < 1L,
    Bin := 1L
  ]

  calibration_dt[
    Bin > number_of_bins,
    Bin := number_of_bins
  ]

  result <- calibration_dt[
    ,
    .(
      n = .N,
      Cases = sum(truth),
      Mean_predicted = mean(probability),
      Observed = mean(truth)
    ),
    by = Bin
  ]

  result[
    ,
    Lower95 := vapply(
      seq_len(.N),
      function(row_index) {
        stats::binom.test(
          x = Cases[[row_index]],
          n = n[[row_index]]
        )$conf.int[[1L]]
      },
      FUN.VALUE = numeric(1)
    )
  ]

  result[
    ,
    Upper95 := vapply(
      seq_len(.N),
      function(row_index) {
        stats::binom.test(
          x = Cases[[row_index]],
          n = n[[row_index]]
        )$conf.int[[2L]]
      },
      FUN.VALUE = numeric(1)
    )
  ]

  result[
    ,
    Group := label
  ]

  result
}


get_pr_precision_on_grid <- function(
  truth,
  probability,
  recall_grid
) {

  truth <- as.character(truth)

  current_pr <- PRROC::pr.curve(
    scores.class0 = probability[
      truth == POSITIVE_CLASS
    ],
    scores.class1 = probability[
      truth == NEGATIVE_CLASS
    ],
    curve = TRUE
  )

  current_curve <- data.table::data.table(
    Recall = as.numeric(
      current_pr$curve[
        ,
        1L
      ]
    ),
    Precision = as.numeric(
      current_pr$curve[
        ,
        2L
      ]
    )
  )

  current_curve <- current_curve[
    is.finite(Recall) &
      is.finite(Precision)
  ]

  current_curve <- current_curve[
    ,
    .(
      Precision = max(
        Precision,
        na.rm = TRUE
      )
    ),
    by = Recall
  ]

  data.table::setorder(
    current_curve,
    Recall
  )

  if (nrow(current_curve) < 2L) {
    return(
      rep(
        NA_real_,
        length(recall_grid)
      )
    )
  }

  stats::approx(
    x = current_curve$Recall,
    y = current_curve$Precision,
    xout = recall_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}


bootstrap_pr_ci <- function(
  truth,
  probability,
  recall_grid = PR_RECALL_GRID,
  bootstraps = PR_BOOTSTRAPS,
  seed = 7832L
) {

  truth <- as.character(truth)
  probability <- as.numeric(probability)

  case_index <- which(
    truth == POSITIVE_CLASS
  )

  control_index <- which(
    truth == NEGATIVE_CLASS
  )

  if (
    length(case_index) < 2L ||
      length(control_index) < 2L
  ) {
    stop(
      "PR bootstrap requires at least two cases and two controls."
    )
  }

  bootstrap_matrix <- matrix(
    NA_real_,
    nrow = bootstraps,
    ncol = length(recall_grid)
  )

  set.seed(seed)

  for (bootstrap_index in seq_len(bootstraps)) {

    sampled_case_index <- sample(
      case_index,
      size = length(case_index),
      replace = TRUE
    )

    sampled_control_index <- sample(
      control_index,
      size = length(control_index),
      replace = TRUE
    )

    sampled_index <- c(
      sampled_case_index,
      sampled_control_index
    )

    bootstrap_matrix[
      bootstrap_index,
    ] <- get_pr_precision_on_grid(
      truth = truth[
        sampled_index
      ],
      probability = probability[
        sampled_index
      ],
      recall_grid = recall_grid
    )
  }

  data.table::data.table(
    Recall = recall_grid,
    Precision_lower = apply(
      bootstrap_matrix,
      2L,
      stats::quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    Precision_median = apply(
      bootstrap_matrix,
      2L,
      stats::quantile,
      probs = 0.500,
      na.rm = TRUE
    ),
    Precision_upper = apply(
      bootstrap_matrix,
      2L,
      stats::quantile,
      probs = 0.975,
      na.rm = TRUE
    )
  )
}

# ======================================================================
# 14. SHARED OUTER FOLDS
# ======================================================================

outer_fold_id <- make_joint_stratified_folds(
  metadata = development_metadata,
  number_of_folds = OUTER_FOLDS,
  seed = SEED_FOUNDATION
)

development_metadata[
  ,
  foundation_outer_fold := outer_fold_id
]

outer_fold_audit <- development_metadata[
  ,
  .(
    n = .N,
    Cases = sum(outcome == POSITIVE_CLASS),
    Controls = sum(outcome == NEGATIVE_CLASS)
  ),
  by = .(
    foundation_outer_fold,
    cohort_group
  )
][
  order(
    foundation_outer_fold,
    cohort_group
  )
]

print(outer_fold_audit)

fwrite(
  outer_fold_audit,
  file.path(
    FOUNDATION_DIR,
    "00_outer_fold_audit.csv"
  )
)


# ======================================================================
# 15. POOLED NESTED CV
# ======================================================================

outer_prediction_rows <- list()
outer_metric_rows <- list()
inner_tuning_rows <- list()
selected_gene_rows <- list()

prediction_counter <- 0L
metric_counter <- 0L
tuning_counter <- 0L
gene_counter <- 0L

for (outer_fold in seq_len(OUTER_FOLDS)) {

  cat(
    "\n\n============================================================\n",
    "CD FOUNDATION EXPERIMENT — OUTER FOLD ",
    outer_fold,
    " OF ",
    OUTER_FOLDS,
    "\n============================================================\n",
    sep = ""
  )

  outer_train_metadata <- development_metadata[
    foundation_outer_fold != outer_fold
  ]

  outer_test_metadata <- development_metadata[
    foundation_outer_fold == outer_fold
  ]

  outer_train_counts <- development_counts[
    ,
    outer_train_metadata$sample_id,
    drop = FALSE
  ]

  outer_test_counts <- development_counts[
    ,
    outer_test_metadata$sample_id,
    drop = FALSE
  ]

  inner_fold_id <- make_joint_stratified_folds(
    metadata = outer_train_metadata,
    number_of_folds = INNER_FOLDS,
    seed = SEED_FOUNDATION +
      outer_fold *
      100L
  )

  outer_train_metadata <- data.table::copy(
    outer_train_metadata
  )

  outer_train_metadata[
    ,
    foundation_inner_fold := inner_fold_id
  ]

  candidate_inner_rows <- list()
  candidate_inner_counter <- 0L

  for (inner_fold in seq_len(INNER_FOLDS)) {

    message(
      "  Inner fold ",
      inner_fold,
      "/",
      INNER_FOLDS
    )

    inner_train_metadata <- outer_train_metadata[
      foundation_inner_fold != inner_fold
    ]

    inner_valid_metadata <- outer_train_metadata[
      foundation_inner_fold == inner_fold
    ]

    inner_train_counts <- development_counts[
      ,
      inner_train_metadata$sample_id,
      drop = FALSE
    ]

    inner_valid_counts <- development_counts[
      ,
      inner_valid_metadata$sample_id,
      drop = FALSE
    ]

    inner_representation <- foundation_fit_apply_logcpm(
      train_counts = inner_train_counts,
      new_counts = inner_valid_counts,
      train_metadata = inner_train_metadata
    )

    inner_genes <- select_foundation_rf_features(
      train_matrix = inner_representation$train,
      train_metadata = inner_train_metadata,
      top_n = FOUNDATION_TOP_N,
      seed = SEED_FOUNDATION +
        outer_fold *
        10000L +
        inner_fold *
        100L
    )

    inner_x_train <- t(
      inner_representation$train[
        inner_genes,
        ,
        drop = FALSE
      ]
    )

    inner_x_valid <- t(
      inner_representation$new[
        inner_genes,
        ,
        drop = FALSE
      ]
    )

    for (config_index in seq_len(nrow(lightgbm_grid))) {

      current_config <- lightgbm_grid[
        config_index
      ]

      current_probability <- fit_predict_lightgbm(
        x_train = inner_x_train,
        y_train = inner_train_metadata$outcome,
        x_test = inner_x_valid,
        parameters = current_config,
        seed = SEED_FOUNDATION +
          outer_fold *
          100000L +
          inner_fold *
          1000L +
          config_index
      )

      current_metric <- evaluate_binary_predictions(
        truth = inner_valid_metadata$outcome,
        probability = current_probability
      )

      candidate_inner_counter <- candidate_inner_counter +
        1L

      candidate_inner_rows[[candidate_inner_counter]] <- data.table::data.table(
        Outer_fold = outer_fold,
        Inner_fold = inner_fold,
        Config_ID = current_config$Config_ID[[1L]],
        AUROC = current_metric$AUROC[[1L]]
      )
    }
  }

  candidate_inner_dt <- data.table::rbindlist(
    candidate_inner_rows
  )

  candidate_summary <- candidate_inner_dt[
    ,
    .(
      Mean_inner_AUROC = mean(AUROC),
      SD_inner_AUROC = stats::sd(AUROC)
    ),
    by = Config_ID
  ]

  data.table::setorder(
    candidate_summary,
    -Mean_inner_AUROC,
    SD_inner_AUROC,
    Config_ID
  )

  best_config_id <- candidate_summary$Config_ID[[1L]]

  best_lightgbm_config <- lightgbm_grid[
    Config_ID == best_config_id
  ]

  cat(
    "Selected LightGBM config: ",
    best_config_id,
    " | mean inner AUROC ",
    sprintf(
      "%.3f",
      candidate_summary$Mean_inner_AUROC[[1L]]
    ),
    "\n",
    sep = ""
  )

  tuning_counter <- tuning_counter +
    1L

  inner_tuning_rows[[tuning_counter]] <- data.table::copy(
    candidate_summary
  )[
    ,
    Outer_fold := outer_fold
  ]

  outer_representation <- foundation_fit_apply_logcpm(
    train_counts = outer_train_counts,
    new_counts = outer_test_counts,
    train_metadata = outer_train_metadata
  )

  outer_genes <- select_foundation_rf_features(
    train_matrix = outer_representation$train,
    train_metadata = outer_train_metadata,
    top_n = FOUNDATION_TOP_N,
    seed = SEED_FOUNDATION +
      outer_fold *
      1000000L
  )

  outer_x_train <- t(
    outer_representation$train[
      outer_genes,
      ,
      drop = FALSE
    ]
  )

  outer_x_test <- t(
    outer_representation$new[
      outer_genes,
      ,
      drop = FALSE
    ]
  )

  gene_counter <- gene_counter +
    1L

  selected_gene_rows[[gene_counter]] <- data.table::data.table(
    Outer_fold = outer_fold,
    Rank = seq_along(outer_genes),
    Gene = outer_genes
  )

  lightgbm_start <- Sys.time()

  lightgbm_probability <- fit_predict_lightgbm(
    x_train = outer_x_train,
    y_train = outer_train_metadata$outcome,
    x_test = outer_x_test,
    parameters = best_lightgbm_config,
    seed = SEED_FOUNDATION +
      outer_fold *
      2000000L
  )

  lightgbm_seconds <- as.numeric(
    difftime(
      Sys.time(),
      lightgbm_start,
      units = "secs"
    )
  )

  tabpfn_start <- Sys.time()

  tabpfn_probability <- fit_predict_tabpfn(
    x_train = outer_x_train,
    y_train = outer_train_metadata$outcome,
    x_test = outer_x_test,
    seed = SEED_FOUNDATION +
      outer_fold *
      3000000L
  )

  tabpfn_seconds <- as.numeric(
    difftime(
      Sys.time(),
      tabpfn_start,
      units = "secs"
    )
  )

  for (model_name in MODEL_ORDER) {

    current_probability <- if (
      model_name == "LightGBM"
    ) {
      lightgbm_probability
    } else {
      tabpfn_probability
    }

    current_seconds <- if (
      model_name == "LightGBM"
    ) {
      lightgbm_seconds
    } else {
      tabpfn_seconds
    }

    prediction_counter <- prediction_counter +
      1L

    outer_prediction_rows[[prediction_counter]] <- data.table::data.table(
      Outer_fold = outer_fold,
      Model = model_name,
      sample_id = outer_test_metadata$sample_id,
      cohort_group = outer_test_metadata$cohort_group,
      truth = as.character(outer_test_metadata$outcome),
      probability_case = current_probability
    )

    current_metrics <- evaluate_binary_predictions(
      truth = outer_test_metadata$outcome,
      probability = current_probability
    )

    metric_counter <- metric_counter +
      1L

    outer_metric_rows[[metric_counter]] <- data.table::data.table(
      Outer_fold = outer_fold,
      Model = model_name,
      AUROC = current_metrics$AUROC[[1L]],
      AUPRC = current_metrics$AUPRC[[1L]],
      BACC = current_metrics$BACC[[1L]],
      Sensitivity = current_metrics$Sensitivity[[1L]],
      Specificity = current_metrics$Specificity[[1L]],
      Brier = current_metrics$Brier[[1L]],
      Runtime_seconds = current_seconds,
      LightGBM_Config = ifelse(
        model_name == "LightGBM",
        best_config_id,
        NA_character_
      )
    )
  }
}


# ======================================================================
# 16. COMBINE AND SAVE NESTED-CV RESULTS
# ======================================================================

foundation_predictions <- data.table::rbindlist(
  outer_prediction_rows,
  use.names = TRUE,
  fill = TRUE
)

foundation_fold_metrics <- data.table::rbindlist(
  outer_metric_rows,
  use.names = TRUE,
  fill = TRUE
)

foundation_inner_tuning <- data.table::rbindlist(
  inner_tuning_rows,
  use.names = TRUE,
  fill = TRUE
)

foundation_selected_genes <- data.table::rbindlist(
  selected_gene_rows,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  foundation_predictions,
  file.path(
    FOUNDATION_DIR,
    "01_nestedCV_predictions.csv"
  )
)

fwrite(
  foundation_fold_metrics,
  file.path(
    FOUNDATION_DIR,
    "02_nestedCV_fold_metrics.csv"
  )
)

fwrite(
  foundation_inner_tuning,
  file.path(
    FOUNDATION_DIR,
    "03_LightGBM_inner_tuning.csv"
  )
)

fwrite(
  foundation_selected_genes,
  file.path(
    FOUNDATION_DIR,
    "04_outer_fold_selected_genes.csv"
  )
)


# ======================================================================
# 17. NESTED-CV SUMMARY
# ======================================================================

foundation_fold_summary <- foundation_fold_metrics[
  ,
  .(
    Outer_folds = .N,
    Mean_AUROC = mean(AUROC),
    SD_AUROC = stats::sd(AUROC),
    Worst_AUROC = min(AUROC),
    Best_AUROC = max(AUROC),
    Mean_AUPRC = mean(AUPRC),
    Mean_BACC = mean(BACC),
    Mean_Brier = mean(Brier),
    Total_runtime_seconds = sum(Runtime_seconds)
  ),
  by = Model
]

data.table::setorder(
  foundation_fold_summary,
  -Mean_AUROC
)

print(foundation_fold_summary)

fwrite(
  foundation_fold_summary,
  file.path(
    FOUNDATION_DIR,
    "05_model_fold_summary.csv"
  )
)


# ======================================================================
# 18. POOLED OUT-OF-FOLD PERFORMANCE
# ======================================================================

pooled_rows <- list()

for (model_index in seq_along(MODEL_ORDER)) {

  current_model <- MODEL_ORDER[[model_index]]

  current_dt <- foundation_predictions[
    Model == current_model
  ]

  stopifnot(
    data.table::uniqueN(current_dt$sample_id) ==
      nrow(development_metadata)
  )

  current_roc <- pROC::roc(
    response = factor(
      current_dt$truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    predictor = current_dt$probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  current_ci <- as.numeric(
    pROC::ci.auc(
      current_roc,
      method = "delong",
      conf.level = 0.95
    )
  )

  current_pr <- PRROC::pr.curve(
    scores.class0 = current_dt$probability_case[
      current_dt$truth == POSITIVE_CLASS
    ],
    scores.class1 = current_dt$probability_case[
      current_dt$truth == NEGATIVE_CLASS
    ],
    curve = FALSE
  )

  pooled_rows[[model_index]] <- data.table::data.table(
    Model = current_model,
    n = nrow(current_dt),
    AUROC = as.numeric(pROC::auc(current_roc)),
    AUROC_lower95 = current_ci[[1L]],
    AUROC_upper95 = current_ci[[3L]],
    AUPRC = current_pr$auc.integral,
    Brier = mean(
      (
        current_dt$probability_case -
          as.integer(
            current_dt$truth == POSITIVE_CLASS
          )
      )^2
    )
  )
}

foundation_pooled_summary <- data.table::rbindlist(
  pooled_rows
)

print(foundation_pooled_summary)

fwrite(
  foundation_pooled_summary,
  file.path(
    FOUNDATION_DIR,
    "06_pooled_OOF_performance.csv"
  )
)


# ======================================================================
# 19. PAIRED OOF DELONG COMPARISON
# ======================================================================

lightgbm_oof <- foundation_predictions[
  Model == "LightGBM"
][
  order(sample_id)
]

tabpfn_oof <- foundation_predictions[
  Model == FOUNDATION_LABEL
][
  order(sample_id)
]

stopifnot(
  identical(
    lightgbm_oof$sample_id,
    tabpfn_oof$sample_id
  )
)

stopifnot(
  identical(
    lightgbm_oof$truth,
    tabpfn_oof$truth
  )
)

lightgbm_roc <- pROC::roc(
  response = factor(
    lightgbm_oof$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),
  predictor = lightgbm_oof$probability_case,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)

tabpfn_roc <- pROC::roc(
  response = factor(
    tabpfn_oof$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),
  predictor = tabpfn_oof$probability_case,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)

paired_delong <- pROC::roc.test(
  lightgbm_roc,
  tabpfn_roc,
  method = "delong",
  paired = TRUE
)

lightgbm_auc <- as.numeric(
  pROC::auc(lightgbm_roc)
)

tabpfn_auc <- as.numeric(
  pROC::auc(tabpfn_roc)
)

foundation_comparison <- data.table::data.table(
  LightGBM_AUROC = lightgbm_auc,
  TabPFN_AUROC = tabpfn_auc,
  Delta_AUROC_TabPFN_minus_LightGBM = tabpfn_auc - lightgbm_auc,
  Paired_DeLong_p = as.numeric(paired_delong$p.value)
)

print(foundation_comparison)

fwrite(
  foundation_comparison,
  file.path(
    FOUNDATION_DIR,
    "07_paired_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 20. PUBLICATION NESTED-CV TABLE
# ======================================================================

publication_foundation_table <- merge(
  foundation_pooled_summary,
  foundation_fold_summary[
    ,
    .(
      Model,
      Mean_outer_AUROC = Mean_AUROC,
      SD_outer_AUROC = SD_AUROC,
      Worst_outer_AUROC = Worst_AUROC,
      Runtime_seconds = Total_runtime_seconds
    )
  ],
  by = "Model"
)

publication_foundation_table[
  ,
  `AUROC (95% CI)` := sprintf(
    "%.3f (%.3f–%.3f)",
    AUROC,
    AUROC_lower95,
    AUROC_upper95
  )
]

publication_foundation_table[
  ,
  `Mean outer AUROC (SD)` := sprintf(
    "%.3f (%.3f)",
    Mean_outer_AUROC,
    SD_outer_AUROC
  )
]

publication_foundation_table[
  ,
  `Worst outer AUROC` := sprintf(
    "%.3f",
    Worst_outer_AUROC
  )
]

publication_foundation_table[
  ,
  AUPRC := round(AUPRC, 3L)
]

publication_foundation_table[
  ,
  Brier := round(Brier, 3L)
]

publication_foundation_table <- publication_foundation_table[
  ,
  .(
    Model,
    n,
    `AUROC (95% CI)`,
    `Mean outer AUROC (SD)`,
    `Worst outer AUROC`,
    AUPRC,
    Brier,
    Runtime_seconds
  )
]

print(publication_foundation_table)

fwrite(
  publication_foundation_table,
  file.path(
    FOUNDATION_DIR,
    "Table_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 21. GENERIC PUBLICATION COMPARISON PANEL
# ======================================================================

make_model_comparison_panel <- function(
  prediction_data,
  panel_title,
  output_dir,
  file_stem,
  seed_offset = 0L
) {

  prediction_data <- data.table::copy(
    prediction_data
  )

  prediction_data[
    ,
    Model := factor(
      as.character(Model),
      levels = MODEL_ORDER
    )
  ]

  roc_curve_rows <- list()
  roc_ribbon_rows <- list()
  roc_labels <- character()

  pr_curve_rows <- list()
  pr_ribbon_rows <- list()
  pr_labels <- character()

  calibration_rows <- list()

  for (model_index in seq_along(MODEL_ORDER)) {

    current_model <- MODEL_ORDER[[model_index]]

    current_dt <- prediction_data[
      Model == current_model
    ]

    current_roc <- pROC::roc(
      response = factor(
        current_dt$truth,
        levels = c(
          NEGATIVE_CLASS,
          POSITIVE_CLASS
        )
      ),
      predictor = current_dt$probability_case,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      ),
      direction = "<",
      quiet = TRUE
    )

    current_auc <- as.numeric(
      pROC::auc(current_roc)
    )

    current_ci <- as.numeric(
      pROC::ci.auc(
        current_roc,
        method = "delong",
        conf.level = 0.95
      )
    )

    current_coordinates <- pROC::coords(
      current_roc,
      x = "all",
      ret = c(
        "specificity",
        "sensitivity"
      ),
      transpose = FALSE
    )

    current_curve_dt <- data.table::data.table(
      Model = current_model,
      FPR = 1 - as.numeric(current_coordinates$specificity),
      TPR = as.numeric(current_coordinates$sensitivity)
    )

    data.table::setorder(
      current_curve_dt,
      FPR,
      TPR
    )

    roc_curve_rows[[model_index]] <- current_curve_dt

    set.seed(
      SEED_FOUNDATION +
        seed_offset +
        model_index *
        10000L
    )

    current_curve_ci <- pROC::ci.se(
      current_roc,
      specificities = ROC_SPECIFICITY_GRID,
      conf.level = 0.95,
      boot.n = MODEL_ROC_BOOTSTRAPS,
      boot.stratified = TRUE
    )

    roc_ribbon_rows[[model_index]] <- data.table::data.table(
      Model = current_model,
      FPR = 1 - ROC_SPECIFICITY_GRID,
      Lower95 = as.numeric(current_curve_ci[, 1L]),
      Upper95 = as.numeric(current_curve_ci[, 3L])
    )

    roc_labels[[current_model]] <- paste0(
      current_model,
      ": ",
      sprintf("%.3f", current_auc),
      " (",
      sprintf("%.3f", current_ci[[1L]]),
      "–",
      sprintf("%.3f", current_ci[[3L]]),
      ")"
    )

    current_pr <- PRROC::pr.curve(
      scores.class0 = current_dt$probability_case[
        current_dt$truth == POSITIVE_CLASS
      ],
      scores.class1 = current_dt$probability_case[
        current_dt$truth == NEGATIVE_CLASS
      ],
      curve = TRUE
    )

    pr_curve_rows[[model_index]] <- data.table::data.table(
      Model = current_model,
      Recall = as.numeric(current_pr$curve[, 1L]),
      Precision = as.numeric(current_pr$curve[, 2L])
    )

    current_pr_ribbon <- bootstrap_pr_ci(
      truth = current_dt$truth,
      probability = current_dt$probability_case,
      recall_grid = PR_RECALL_GRID,
      bootstraps = PR_BOOTSTRAPS,
      seed = SEED_FOUNDATION +
        seed_offset +
        model_index *
        20000L
    )

    current_pr_ribbon[
      ,
      Model := current_model
    ]

    pr_ribbon_rows[[model_index]] <- current_pr_ribbon

    pr_labels[[current_model]] <- paste0(
      current_model,
      ": AUPRC ",
      sprintf(
        "%.3f",
        current_pr$auc.integral
      )
    )

    calibration_rows[[model_index]] <- make_calibration_data(
      truth = current_dt$truth,
      probability = current_dt$probability_case,
      label = current_model,
      number_of_bins = 8L
    )
  }

  roc_curve_dt <- data.table::rbindlist(
    roc_curve_rows
  )

  roc_ribbon_dt <- data.table::rbindlist(
    roc_ribbon_rows
  )

  pr_curve_dt <- data.table::rbindlist(
    pr_curve_rows
  )

  pr_ribbon_dt <- data.table::rbindlist(
    pr_ribbon_rows
  )

  calibration_dt <- data.table::rbindlist(
    calibration_rows
  )

  for (current_dt_name in c(
    "roc_curve_dt",
    "roc_ribbon_dt",
    "pr_curve_dt",
    "pr_ribbon_dt"
  )) {
    current_dt <- get(current_dt_name)
    current_dt[
      ,
      Model := factor(
        Model,
        levels = MODEL_ORDER
      )
    ]
    assign(current_dt_name, current_dt)
  }

  calibration_dt[
    ,
    Group := factor(
      Group,
      levels = MODEL_ORDER
    )
  ]

  paired_dt <- data.table::dcast(
    prediction_data,
    sample_id + truth ~ Model,
    value.var = "probability_case"
  )

  paired_truth <- factor(
    paired_dt$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  paired_lgbm_roc <- pROC::roc(
    response = paired_truth,
    predictor = paired_dt$LightGBM,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  paired_tabpfn_roc <- pROC::roc(
    response = paired_truth,
    predictor = paired_dt[[FOUNDATION_LABEL]],
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  paired_delong_local <- pROC::roc.test(
    paired_lgbm_roc,
    paired_tabpfn_roc,
    method = "delong",
    paired = TRUE
  )

  delta_auc_local <- as.numeric(
    pROC::auc(paired_tabpfn_roc) -
      pROC::auc(paired_lgbm_roc)
  )

  delong_p_local <- as.numeric(
    paired_delong_local$p.value
  )

  comparison_annotation <- paste0(
    "ΔAUROC (TabPFN - LightGBM) = ",
    sprintf("%+.3f", delta_auc_local),
    "\n",
    ifelse(
      delong_p_local < 0.001,
      "Paired DeLong p < 0.001",
      paste0(
        "Paired DeLong p = ",
        sprintf("%.3f", delong_p_local)
      )
    )
  )

  comparison_prevalence <- mean(
    prediction_data$truth == POSITIVE_CLASS
  )

  p_model_roc <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = roc_ribbon_dt,
      ggplot2::aes(
        x = FPR,
        ymin = Lower95,
        ymax = Upper95,
        fill = Model,
        group = Model
      ),
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.65
    ) +
    ggplot2::geom_step(
      data = roc_curve_dt,
      ggplot2::aes(
        x = FPR,
        y = TPR,
        colour = Model,
        group = Model
      ),
      direction = "vh",
      linewidth = 1.15
    ) +
    ggplot2::annotate(
      geom = "label",
      x = 0.97,
      y = 0.05,
      label = comparison_annotation,
      hjust = 1,
      vjust = 0,
      size = 3.2,
      fill = "white"
    ) +
    ggplot2::scale_colour_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      labels = roc_labels[MODEL_ORDER],
      name = "Model — AUROC (95% CI)"
    ) +
    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Discrimination",
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    publication_theme() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left"
    )

  p_model_pr <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pr_ribbon_dt,
      ggplot2::aes(
        x = Recall,
        ymin = Precision_lower,
        ymax = Precision_upper,
        fill = Model,
        group = Model
      ),
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_hline(
      yintercept = comparison_prevalence,
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.70
    ) +
    ggplot2::geom_path(
      data = pr_curve_dt,
      ggplot2::aes(
        x = Recall,
        y = Precision,
        colour = Model,
        group = Model
      ),
      linewidth = 1.10
    ) +
    ggplot2::scale_colour_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      labels = pr_labels[MODEL_ORDER],
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Precision-recall",
      subtitle = paste0(
        "Dashed line = prevalence ",
        sprintf("%.3f", comparison_prevalence)
      ),
      x = "Recall",
      y = "Precision"
    ) +
    publication_theme() +
    ggplot2::theme(
      legend.position = "top"
    )

  p_model_calibration <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = calibration_dt,
      ggplot2::aes(
        x = Mean_predicted,
        ymin = Lower95,
        ymax = Upper95,
        fill = Group,
        group = Group
      ),
      alpha = 0.13,
      colour = NA
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::geom_line(
      data = calibration_dt,
      ggplot2::aes(
        x = Mean_predicted,
        y = Observed,
        colour = Group,
        group = Group
      ),
      linewidth = 0.95
    ) +
    ggplot2::geom_point(
      data = calibration_dt,
      ggplot2::aes(
        x = Mean_predicted,
        y = Observed,
        colour = Group
      ),
      size = 2.5
    ) +
    ggplot2::scale_colour_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      name = "Model"
    ) +
    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      guide = "none"
    ) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Calibration",
      x = "Mean predicted probability",
      y = "Observed case proportion"
    ) +
    publication_theme() +
    ggplot2::theme(
      legend.position = "top"
    )

  p_main <- (
    p_model_roc +
      p_model_pr +
      p_model_calibration
  ) +
    patchwork::plot_layout(
      ncol = 3
    ) +
    patchwork::plot_annotation(
      title = panel_title,
      subtitle = "Identical samples, representation and training-selected transcriptomic features",
      tag_levels = "A"
    )

  print(p_main)

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(file_stem, ".pdf")
    ),
    plot = p_main,
    width = 15.8,
    height = 5.4,
    units = "in",
    dpi = 600,
    device = grDevices::cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(file_stem, ".png")
    ),
    plot = p_main,
    width = 15.8,
    height = 5.4,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  fwrite(
    roc_curve_dt,
    file.path(
      output_dir,
      paste0(file_stem, "_ROC_coordinates.csv")
    )
  )

  fwrite(
    roc_ribbon_dt,
    file.path(
      output_dir,
      paste0(file_stem, "_ROC_95CI.csv")
    )
  )

  fwrite(
    pr_ribbon_dt,
    file.path(
      output_dir,
      paste0(file_stem, "_PR_95CI.csv")
    )
  )

  fwrite(
    calibration_dt,
    file.path(
      output_dir,
      paste0(file_stem, "_calibration.csv")
    )
  )

  list(
    plot = p_main,
    paired_delong = paired_delong_local,
    delta_auc = delta_auc_local,
    p_value = delong_p_local
  )
}


# ======================================================================
# 22. NESTED-CV PUBLICATION FIGURES
# ======================================================================

nested_comparison_panel <- make_model_comparison_panel(
  prediction_data = foundation_predictions,
  panel_title = "Task-specific gradient boosting versus a tabular foundation model in CD",
  output_dir = PUBLICATION_DIR,
  file_stem = "Figure_MAIN_CD_LightGBM_vs_TabPFN_nestedCV_ROC_PR_calibration",
  seed_offset = 0L
)

paired_fold_dt <- data.table::copy(
  foundation_fold_metrics
)

paired_fold_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

p_paired_outer_auc <- ggplot2::ggplot(
  paired_fold_dt,
  ggplot2::aes(
    x = Model,
    y = AUROC,
    group = Outer_fold
  )
) +
  ggplot2::geom_line(
    colour = "grey65",
    linewidth = 0.70
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = Model
    ),
    size = 3.2
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey45"
  ) +
  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0.45, 1)
  ) +
  ggplot2::labs(
    title = "Paired outer-fold discrimination in CD",
    subtitle = "Lines connect predictions from the identical nested-CV outer fold",
    x = NULL,
    y = "Outer-fold AUROC"
  ) +
  publication_theme()

print(p_paired_outer_auc)

ggplot2::ggsave(
  filename = file.path(
    PUBLICATION_DIR,
    "Figure_S1_CD_paired_outer_fold_AUROC.pdf"
  ),
  plot = p_paired_outer_auc,
  width = 6.4,
  height = 5.2,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


# ======================================================================
# 23. STRICT LOCO: SOURCE-ONLY LIGHTGBM TUNING
# ======================================================================

LOCO_TOP_N <- FOUNDATION_TOP_N
LOCO_INNER_FOLDS <- 3L


tune_lightgbm_source_only <- function(
  source_counts,
  source_metadata,
  target_index
) {

  inner_fold_id <- make_joint_stratified_folds(
    metadata = source_metadata,
    number_of_folds = LOCO_INNER_FOLDS,
    seed = SEED_FOUNDATION +
      target_index *
      100000L
  )

  source_metadata <- data.table::copy(
    source_metadata
  )

  source_metadata[
    ,
    loco_inner_fold := inner_fold_id
  ]

  tuning_rows <- list()
  tuning_counter <- 0L

  for (inner_fold in seq_len(LOCO_INNER_FOLDS)) {

    inner_train_metadata <- source_metadata[
      loco_inner_fold != inner_fold
    ]

    inner_valid_metadata <- source_metadata[
      loco_inner_fold == inner_fold
    ]

    inner_train_counts <- source_counts[
      ,
      inner_train_metadata$sample_id,
      drop = FALSE
    ]

    inner_valid_counts <- source_counts[
      ,
      inner_valid_metadata$sample_id,
      drop = FALSE
    ]

    inner_representation <- foundation_fit_apply_logcpm(
      train_counts = inner_train_counts,
      new_counts = inner_valid_counts,
      train_metadata = inner_train_metadata
    )

    inner_genes <- select_foundation_rf_features(
      train_matrix = inner_representation$train,
      train_metadata = inner_train_metadata,
      top_n = LOCO_TOP_N,
      seed = SEED_FOUNDATION +
        target_index *
        1000000L +
        inner_fold *
        1000L
    )

    inner_x_train <- t(
      inner_representation$train[
        inner_genes,
        ,
        drop = FALSE
      ]
    )

    inner_x_valid <- t(
      inner_representation$new[
        inner_genes,
        ,
        drop = FALSE
      ]
    )

    for (config_index in seq_len(nrow(lightgbm_grid))) {

      current_config <- lightgbm_grid[
        config_index
      ]

      probability <- fit_predict_lightgbm(
        x_train = inner_x_train,
        y_train = inner_train_metadata$outcome,
        x_test = inner_x_valid,
        parameters = current_config,
        seed = SEED_FOUNDATION +
          target_index *
          10000000L +
          inner_fold *
          10000L +
          config_index
      )

      metric <- evaluate_binary_predictions(
        truth = inner_valid_metadata$outcome,
        probability = probability
      )

      tuning_counter <- tuning_counter +
        1L

      tuning_rows[[tuning_counter]] <- data.table::data.table(
        Inner_fold = inner_fold,
        Config_ID = current_config$Config_ID[[1L]],
        AUROC = metric$AUROC[[1L]]
      )
    }
  }

  tuning_dt <- data.table::rbindlist(
    tuning_rows
  )

  tuning_summary <- tuning_dt[
    ,
    .(
      Mean_inner_AUROC = mean(AUROC),
      SD_inner_AUROC = stats::sd(AUROC)
    ),
    by = Config_ID
  ]

  data.table::setorder(
    tuning_summary,
    -Mean_inner_AUROC,
    SD_inner_AUROC,
    Config_ID
  )

  best_config_id <- tuning_summary$Config_ID[[1L]]

  list(
    best_config = lightgbm_grid[
      Config_ID == best_config_id
    ],
    summary = tuning_summary,
    raw = tuning_dt
  )
}


# ======================================================================
# 24. STRICT LOCO LOOP
# ======================================================================

loco_prediction_rows <- list()
loco_metric_rows <- list()
loco_tuning_rows <- list()
loco_gene_rows <- list()

prediction_counter <- 0L
metric_counter <- 0L
tuning_counter <- 0L
gene_counter <- 0L

for (target_index in seq_along(COHORT_ORDER)) {

  target_cohort <- COHORT_ORDER[[target_index]]

  cat(
    "\n\n============================================================\n",
    "CD STRICT LOCO: ",
    target_cohort,
    " COMPLETELY HELD OUT\n",
    "============================================================\n",
    sep = ""
  )

  source_metadata <- development_metadata[
    cohort_group != target_cohort
  ]

  target_metadata <- development_metadata[
    cohort_group == target_cohort
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

  lightgbm_tuning <- tune_lightgbm_source_only(
    source_counts = source_counts,
    source_metadata = source_metadata,
    target_index = target_index
  )

  best_lgbm_config <- lightgbm_tuning$best_config

  cat(
    "Selected LightGBM config: ",
    best_lgbm_config$Config_ID[[1L]],
    "\n",
    sep = ""
  )

  tuning_counter <- tuning_counter +
    1L

  loco_tuning_rows[[tuning_counter]] <- data.table::copy(
    lightgbm_tuning$summary
  )[
    ,
    Target_cohort := target_cohort
  ]

  loco_representation <- foundation_fit_apply_logcpm(
    train_counts = source_counts,
    new_counts = target_counts,
    train_metadata = source_metadata
  )

  loco_genes <- select_foundation_rf_features(
    train_matrix = loco_representation$train,
    train_metadata = source_metadata,
    top_n = LOCO_TOP_N,
    seed = SEED_FOUNDATION +
      target_index *
      20000000L
  )

  loco_x_train <- t(
    loco_representation$train[
      loco_genes,
      ,
      drop = FALSE
    ]
  )

  loco_x_test <- t(
    loco_representation$new[
      loco_genes,
      ,
      drop = FALSE
    ]
  )

  gene_counter <- gene_counter +
    1L

  loco_gene_rows[[gene_counter]] <- data.table::data.table(
    Target_cohort = target_cohort,
    Rank = seq_along(loco_genes),
    Gene = loco_genes
  )

  lgbm_probability <- fit_predict_lightgbm(
    x_train = loco_x_train,
    y_train = source_metadata$outcome,
    x_test = loco_x_test,
    parameters = best_lgbm_config,
    seed = SEED_FOUNDATION +
      target_index *
      30000000L
  )

  tabpfn_probability <- fit_predict_tabpfn(
    x_train = loco_x_train,
    y_train = source_metadata$outcome,
    x_test = loco_x_test,
    seed = SEED_FOUNDATION +
      target_index *
      40000000L
  )

  for (current_model in MODEL_ORDER) {

    probability <- if (
      current_model == "LightGBM"
    ) {
      lgbm_probability
    } else {
      tabpfn_probability
    }

    current_roc <- pROC::roc(
      response = factor(
        target_metadata$outcome,
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

    current_ci <- as.numeric(
      pROC::ci.auc(
        current_roc,
        method = "delong",
        conf.level = 0.95
      )
    )

    current_metrics <- evaluate_binary_predictions(
      truth = target_metadata$outcome,
      probability = probability
    )

    prediction_counter <- prediction_counter +
      1L

    loco_prediction_rows[[prediction_counter]] <- data.table::data.table(
      Target_cohort = target_cohort,
      Model = current_model,
      sample_id = target_metadata$sample_id,
      truth = as.character(target_metadata$outcome),
      probability_case = probability
    )

    metric_counter <- metric_counter +
      1L

    loco_metric_rows[[metric_counter]] <- data.table::data.table(
      Target_cohort = target_cohort,
      Model = current_model,
      n = nrow(target_metadata),
      Cases = sum(
        target_metadata$outcome == POSITIVE_CLASS
      ),
      Controls = sum(
        target_metadata$outcome == NEGATIVE_CLASS
      ),
      AUROC = as.numeric(pROC::auc(current_roc)),
      AUROC_lower95 = current_ci[[1L]],
      AUROC_upper95 = current_ci[[3L]],
      AUPRC = current_metrics$AUPRC[[1L]],
      BACC = current_metrics$BACC[[1L]],
      Sensitivity = current_metrics$Sensitivity[[1L]],
      Specificity = current_metrics$Specificity[[1L]],
      Brier = current_metrics$Brier[[1L]]
    )
  }
}

loco_foundation_predictions <- data.table::rbindlist(
  loco_prediction_rows
)

loco_foundation_metrics <- data.table::rbindlist(
  loco_metric_rows
)

loco_foundation_tuning <- data.table::rbindlist(
  loco_tuning_rows
)

loco_foundation_genes <- data.table::rbindlist(
  loco_gene_rows
)


generalisability_summary <- loco_foundation_metrics[
  ,
  .(
    Mean_LOCO_AUROC = mean(AUROC),
    Worst_cohort_AUROC = min(AUROC),
    SD_LOCO_AUROC = stats::sd(AUROC),
    Mean_AUPRC = mean(AUPRC),
    Mean_BACC = mean(BACC),
    Mean_Brier = mean(Brier)
  ),
  by = Model
]

data.table::setorder(
  generalisability_summary,
  -Worst_cohort_AUROC,
  -Mean_LOCO_AUROC
)

generalisability_summary[
  ,
  Generalisability_rank := seq_len(.N)
]

nested_reference <- foundation_pooled_summary[
  ,
  .(
    Model,
    Nested_OOF_AUROC = AUROC
  )
]

transportability_table <- merge(
  generalisability_summary,
  nested_reference,
  by = "Model"
)

transportability_table[
  ,
  Delta_LOCO_minus_nested :=
    Mean_LOCO_AUROC -
    Nested_OOF_AUROC
]

transportability_table[
  ,
  AUROC_retention_percent :=
    100 *
    Mean_LOCO_AUROC /
    Nested_OOF_AUROC
]

data.table::setorder(
  transportability_table,
  -Worst_cohort_AUROC
)

print(loco_foundation_metrics)
print(generalisability_summary)
print(transportability_table)

fwrite(
  loco_foundation_predictions,
  file.path(
    GENERALISATION_DIR,
    "01_strict_LOCO_predictions.csv"
  )
)

fwrite(
  loco_foundation_metrics,
  file.path(
    GENERALISATION_DIR,
    "02_strict_LOCO_by_cohort_metrics.csv"
  )
)

fwrite(
  generalisability_summary,
  file.path(
    GENERALISATION_DIR,
    "03_generalisability_summary.csv"
  )
)

fwrite(
  transportability_table,
  file.path(
    GENERALISATION_DIR,
    "04_nestedCV_to_LOCO_retention.csv"
  )
)

fwrite(
  loco_foundation_tuning,
  file.path(
    GENERALISATION_DIR,
    "05_source_only_LightGBM_tuning.csv"
  )
)

fwrite(
  loco_foundation_genes,
  file.path(
    GENERALISATION_DIR,
    "06_LOCO_selected_genes.csv"
  )
)


# ======================================================================
# 25. STRICT LOCO TRANSPORTABILITY FIGURE
# ======================================================================

loco_plot_dt <- data.table::copy(
  loco_foundation_metrics
)

loco_plot_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = COHORT_ORDER
  )
]

loco_plot_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

p_generalisation <- ggplot2::ggplot(
  loco_plot_dt,
  ggplot2::aes(
    x = Target_cohort,
    y = AUROC,
    colour = Model,
    group = Model
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey45"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = AUROC_lower95,
      ymax = AUROC_upper95
    ),
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    width = 0.12,
    linewidth = 0.60
  ) +
  ggplot2::geom_line(
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    linewidth = 0.80
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    size = 3.2
  ) +
  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    name = "Model"
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0.45, 1)
  ) +
  ggplot2::labs(
    title = "Transportability to completely unseen CD cohorts",
    subtitle = "Each target cohort was excluded from preprocessing, feature selection, LightGBM tuning and model fitting",
    x = "Completely held-out cohort",
    y = "AUROC (95% CI)"
  ) +
  publication_theme() +
  ggplot2::theme(
    legend.position = "top"
  )

print(p_generalisation)

ggplot2::ggsave(
  filename = file.path(
    GENERALISATION_DIR,
    "Figure_STRICT_LOCO_CD_LightGBM_vs_TabPFN.pdf"
  ),
  plot = p_generalisation,
  width = 7.8,
  height = 5.5,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    GENERALISATION_DIR,
    "Figure_STRICT_LOCO_CD_LightGBM_vs_TabPFN.png"
  ),
  plot = p_generalisation,
  width = 7.8,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 26. SELECT FINAL LIGHTGBM CONFIGURATION FROM INNER-CV ONLY
# ======================================================================

final_lgbm_tuning_summary <- foundation_inner_tuning[
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

data.table::setorder(
  final_lgbm_tuning_summary,
  -Mean_inner_AUROC,
  SD_inner_AUROC,
  -Worst_outer_inner_AUROC,
  Config_ID
)

FINAL_LGBM_CONFIG_ID <- final_lgbm_tuning_summary$Config_ID[[1L]]

FINAL_LGBM_CONFIG <- lightgbm_grid[
  Config_ID == FINAL_LGBM_CONFIG_ID
]

print(final_lgbm_tuning_summary)

cat(
  "\nFinal LightGBM config: ",
  FINAL_LGBM_CONFIG_ID,
  "\n",
  sep = ""
)

fwrite(
  final_lgbm_tuning_summary,
  file.path(
    ICE_COMPARISON_DIR,
    "01_final_LightGBM_config_selection.csv"
  )
)


# ======================================================================
# 27. FINAL DEVELOPMENT-ONLY REPRESENTATION + FEATURE SELECTION
# ======================================================================

ice_representation <- foundation_fit_apply_logcpm(
  train_counts = development_counts,
  new_counts = final_ice_counts,
  train_metadata = development_metadata
)

ICE_SELECTED_GENES <- select_foundation_rf_features(
  train_matrix = ice_representation$train,
  train_metadata = development_metadata,
  top_n = FOUNDATION_TOP_N,
  seed = SEED_FOUNDATION +
    99000001L
)

fwrite(
  data.table::data.table(
    Rank = seq_along(ICE_SELECTED_GENES),
    Gene = ICE_SELECTED_GENES
  ),
  file.path(
    ICE_COMPARISON_DIR,
    "02_final_100_selected_genes.csv"
  )
)

ice_x_train <- t(
  ice_representation$train[
    ICE_SELECTED_GENES,
    ,
    drop = FALSE
  ]
)

ice_x_test <- t(
  ice_representation$new[
    ICE_SELECTED_GENES,
    ,
    drop = FALSE
  ]
)

stopifnot(
  nrow(ice_x_train) ==
    nrow(development_metadata)
)

stopifnot(
  nrow(ice_x_test) ==
    nrow(final_ice_metadata)
)

stopifnot(
  identical(
    rownames(ice_x_test),
    final_ice_metadata$sample_id
  )
)

saveRDS(
  ice_x_train,
  file.path(
    ICE_COMPARISON_DIR,
    "CD_foundation_final_x_train.rds"
  )
)

saveRDS(
  ice_x_test,
  file.path(
    ICE_COMPARISON_DIR,
    "CD_foundation_final_x_ice.rds"
  )
)


# ======================================================================
# 28. FIT FINAL LIGHTGBM AND TABPFN OBJECTS
# ======================================================================

cat("\nFitting final CD LightGBM...\n")

lgbm_start_time <- Sys.time()

FINAL_LGBM_MODEL <- fit_lightgbm_object(
  x_train = ice_x_train,
  y_train = development_metadata$outcome,
  parameters = FINAL_LGBM_CONFIG,
  seed = SEED_FOUNDATION +
    99000002L
)

final_ice_lgbm_probability <- predict_python_classifier(
  model = FINAL_LGBM_MODEL,
  x_test = ice_x_test
)

lgbm_runtime_seconds <- as.numeric(
  difftime(
    Sys.time(),
    lgbm_start_time,
    units = "secs"
  )
)

cat(
  "\nFitting final CD ",
  FOUNDATION_LABEL,
  "...\n",
  sep = ""
)

tabpfn_start_time <- Sys.time()

FINAL_TABPFN_MODEL <- fit_tabpfn_object(
  x_train = ice_x_train,
  y_train = development_metadata$outcome,
  seed = SEED_FOUNDATION +
    99000003L
)

final_ice_tabpfn_probability <- predict_python_classifier(
  model = FINAL_TABPFN_MODEL,
  x_test = ice_x_test
)

tabpfn_runtime_seconds <- as.numeric(
  difftime(
    Sys.time(),
    tabpfn_start_time,
    units = "secs"
  )
)


# ======================================================================
# 29. SAVE FINAL PYTHON MODELS CORRECTLY
# ======================================================================

reticulate::py_save_object(
  FINAL_LGBM_MODEL,
  file.path(
    MODEL_DIR,
    "CD_FINAL_LightGBM_model.pickle"
  )
)

reticulate::py_save_object(
  FINAL_TABPFN_MODEL,
  file.path(
    MODEL_DIR,
    "CD_FINAL_TabPFN_model.pickle"
  )
)

foundation_state <- list(
  disease = "CD",
  representation = "log2(CPM + 1)",
  feature_selection = "RF permutation importance",
  selected_genes = ICE_SELECTED_GENES,
  top_n = FOUNDATION_TOP_N,
  final_lightgbm_config_id = FINAL_LGBM_CONFIG_ID,
  final_lightgbm_config = FINAL_LGBM_CONFIG,
  positive_class = POSITIVE_CLASS,
  negative_class = NEGATIVE_CLASS,
  development_sample_ids = development_metadata$sample_id,
  ice_sample_ids = final_ice_metadata$sample_id,
  cohort_order = COHORT_ORDER,
  tabpfn_device = TABPFN_DEVICE
)

saveRDS(
  foundation_state,
  file.path(
    MODEL_DIR,
    "CD_LightGBM_TabPFN_experiment_state.rds"
  )
)


# ======================================================================
# 30. DIAGNOSTIC FINAL ICE PREDICTIONS
# ======================================================================

ice_model_predictions <- data.table::rbindlist(
  list(
    data.table::data.table(
      Model = "LightGBM",
      sample_id = final_ice_metadata$sample_id,
      cohort_group = final_ice_metadata$cohort_group,
      truth = as.character(final_ice_metadata$outcome),
      probability_case = final_ice_lgbm_probability
    ),
    data.table::data.table(
      Model = FOUNDATION_LABEL,
      sample_id = final_ice_metadata$sample_id,
      cohort_group = final_ice_metadata$cohort_group,
      truth = as.character(final_ice_metadata$outcome),
      probability_case = final_ice_tabpfn_probability
    )
  )
)

ice_model_predictions[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

ice_prediction_audit <- ice_model_predictions[
  ,
  .(
    Rows = .N,
    Unique_samples = data.table::uniqueN(sample_id),
    Cases = sum(truth == POSITIVE_CLASS),
    Controls = sum(truth == NEGATIVE_CLASS)
  ),
  by = Model
]

print(ice_prediction_audit)

fwrite(
  ice_model_predictions,
  file.path(
    ICE_COMPARISON_DIR,
    "03_FINAL_ICE_predictions_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 31. DIAGNOSTIC FINAL ICE PERFORMANCE
# ======================================================================

ICE_THRESHOLD <- 0.50
ice_metric_rows <- list()

for (model_index in seq_along(MODEL_ORDER)) {

  current_model <- MODEL_ORDER[[model_index]]

  current_dt <- ice_model_predictions[
    Model == current_model
  ]

  current_roc <- pROC::roc(
    response = factor(
      current_dt$truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    predictor = current_dt$probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  current_ci <- as.numeric(
    pROC::ci.auc(
      current_roc,
      method = "delong",
      conf.level = 0.95
    )
  )

  current_metrics <- evaluate_binary_predictions(
    truth = current_dt$truth,
    probability = current_dt$probability_case,
    threshold = ICE_THRESHOLD
  )

  ice_metric_rows[[model_index]] <- data.table::data.table(
    Model = current_model,
    n = nrow(current_dt),
    Cases = sum(current_dt$truth == POSITIVE_CLASS),
    Controls = sum(current_dt$truth == NEGATIVE_CLASS),
    AUROC = as.numeric(pROC::auc(current_roc)),
    AUROC_lower95 = current_ci[[1L]],
    AUROC_upper95 = current_ci[[3L]],
    AUPRC = current_metrics$AUPRC[[1L]],
    BACC = current_metrics$BACC[[1L]],
    Sensitivity = current_metrics$Sensitivity[[1L]],
    Specificity = current_metrics$Specificity[[1L]],
    Brier = current_metrics$Brier[[1L]],
    Threshold = ICE_THRESHOLD,
    Runtime_seconds = if (
      current_model == "LightGBM"
    ) {
      lgbm_runtime_seconds
    } else {
      tabpfn_runtime_seconds
    }
  )
}

ice_model_summary <- data.table::rbindlist(
  ice_metric_rows
)

print(ice_model_summary)

fwrite(
  ice_model_summary,
  file.path(
    ICE_COMPARISON_DIR,
    "04_FINAL_ICE_performance_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 32. PAIRED DIAGNOSTIC ICE DELONG
# ======================================================================

paired_ice_dt <- data.table::dcast(
  ice_model_predictions,
  sample_id + cohort_group + truth ~ Model,
  value.var = "probability_case"
)

paired_ice_truth <- factor(
  paired_ice_dt$truth,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  )
)

ice_lgbm_roc <- pROC::roc(
  response = paired_ice_truth,
  predictor = paired_ice_dt$LightGBM,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)

ice_tabpfn_roc <- pROC::roc(
  response = paired_ice_truth,
  predictor = paired_ice_dt[[FOUNDATION_LABEL]],
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)

ice_paired_delong <- pROC::roc.test(
  ice_lgbm_roc,
  ice_tabpfn_roc,
  method = "delong",
  paired = TRUE
)

ice_lgbm_auc <- as.numeric(
  pROC::auc(ice_lgbm_roc)
)

ice_tabpfn_auc <- as.numeric(
  pROC::auc(ice_tabpfn_roc)
)

ice_paired_comparison <- data.table::data.table(
  LightGBM_AUROC = ice_lgbm_auc,
  TabPFN_AUROC = ice_tabpfn_auc,
  Delta_AUROC_TabPFN_minus_LightGBM = ice_tabpfn_auc - ice_lgbm_auc,
  Paired_DeLong_p = as.numeric(ice_paired_delong$p.value)
)

print(ice_paired_comparison)

fwrite(
  ice_paired_comparison,
  file.path(
    ICE_COMPARISON_DIR,
    "05_FINAL_ICE_paired_DeLong_comparison.csv"
  )
)


# ======================================================================
# 33. FINAL ICE PUBLICATION TABLE
# ======================================================================

ice_publication_table <- data.table::copy(
  ice_model_summary
)

ice_publication_table[
  ,
  `AUROC (95% CI)` := sprintf(
    "%.3f (%.3f–%.3f)",
    AUROC,
    AUROC_lower95,
    AUROC_upper95
  )
]

ice_publication_table[
  ,
  `:=`(
    AUPRC = round(AUPRC, 3L),
    BACC = round(BACC, 3L),
    Sensitivity = round(Sensitivity, 3L),
    Specificity = round(Specificity, 3L),
    Brier = round(Brier, 3L),
    Runtime_seconds = round(Runtime_seconds, 2L)
  )
]

ice_publication_table <- ice_publication_table[
  ,
  .(
    Model,
    n,
    Cases,
    Controls,
    `AUROC (95% CI)`,
    AUPRC,
    BACC,
    Sensitivity,
    Specificity,
    Brier,
    Runtime_seconds
  )
]

print(ice_publication_table)

fwrite(
  ice_publication_table,
  file.path(
    ICE_COMPARISON_DIR,
    "Table_FINAL_ICE_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 34. FINAL ICE PUBLICATION ROC + PR + CALIBRATION PANEL
# ======================================================================

ice_comparison_panel <- make_model_comparison_panel(
  prediction_data = ice_model_predictions,
  panel_title = "Diagnostic held-out comparison of LightGBM and TabPFN in CD",
  output_dir = ICE_COMPARISON_DIR,
  file_stem = "Figure_FINAL_ICE_CD_LightGBM_vs_TabPFN_ROC_PR_calibration",
  seed_offset = 500000L
)


# ======================================================================
# 35. PUBLICATION MANIFEST
# ======================================================================

publication_manifest <- data.table::data.table(
  Item = c(
    "Nested-CV predictions",
    "Nested-CV fold metrics",
    "Nested-CV pooled performance",
    "Nested-CV paired DeLong comparison",
    "Nested-CV main ROC/PR/calibration figure",
    "Paired outer-fold AUROC figure",
    "Strict LOCO predictions",
    "Strict LOCO cohort metrics",
    "Strict LOCO summary",
    "Nested-CV to LOCO retention",
    "Strict LOCO transportability figure",
    "Diagnostic ICE predictions",
    "Diagnostic ICE performance",
    "Diagnostic ICE paired DeLong comparison",
    "Diagnostic ICE ROC/PR/calibration figure",
    "Restart-safe LightGBM model",
    "Restart-safe TabPFN model",
    "Experiment state"
  ),
  Path = c(
    file.path(FOUNDATION_DIR, "01_nestedCV_predictions.csv"),
    file.path(FOUNDATION_DIR, "02_nestedCV_fold_metrics.csv"),
    file.path(FOUNDATION_DIR, "06_pooled_OOF_performance.csv"),
    file.path(FOUNDATION_DIR, "07_paired_LightGBM_vs_TabPFN.csv"),
    file.path(PUBLICATION_DIR, "Figure_MAIN_CD_LightGBM_vs_TabPFN_nestedCV_ROC_PR_calibration.pdf"),
    file.path(PUBLICATION_DIR, "Figure_S1_CD_paired_outer_fold_AUROC.pdf"),
    file.path(GENERALISATION_DIR, "01_strict_LOCO_predictions.csv"),
    file.path(GENERALISATION_DIR, "02_strict_LOCO_by_cohort_metrics.csv"),
    file.path(GENERALISATION_DIR, "03_generalisability_summary.csv"),
    file.path(GENERALISATION_DIR, "04_nestedCV_to_LOCO_retention.csv"),
    file.path(GENERALISATION_DIR, "Figure_STRICT_LOCO_CD_LightGBM_vs_TabPFN.pdf"),
    file.path(ICE_COMPARISON_DIR, "03_FINAL_ICE_predictions_LightGBM_vs_TabPFN.csv"),
    file.path(ICE_COMPARISON_DIR, "04_FINAL_ICE_performance_LightGBM_vs_TabPFN.csv"),
    file.path(ICE_COMPARISON_DIR, "05_FINAL_ICE_paired_DeLong_comparison.csv"),
    file.path(ICE_COMPARISON_DIR, "Figure_FINAL_ICE_CD_LightGBM_vs_TabPFN_ROC_PR_calibration.pdf"),
    file.path(MODEL_DIR, "CD_FINAL_LightGBM_model.pickle"),
    file.path(MODEL_DIR, "CD_FINAL_TabPFN_model.pickle"),
    file.path(MODEL_DIR, "CD_LightGBM_TabPFN_experiment_state.rds")
  )
)

fwrite(
  publication_manifest,
  file.path(
    FOUNDATION_DIR,
    "PUBLICATION_MANIFEST.csv"
  )
)


# ======================================================================
# 36. FINAL RESULT AUDIT
# ======================================================================

lgbm_loco <- generalisability_summary[
  Model == "LightGBM"
]

tabpfn_loco <- generalisability_summary[
  Model == FOUNDATION_LABEL
]

cat(
  "\n\n============================================================\n",
  "CD LIGHTGBM vs TABPFN EXPERIMENT COMPLETE\n",
  "============================================================\n",
  "Development n: ", nrow(development_metadata), "\n",
  "Diagnostic ICE n: ", nrow(final_ice_metadata), "\n",
  "Representation: log2(CPM + 1)\n",
  "Shared feature set size: ", FOUNDATION_TOP_N, "\n",
  "TabPFN device: ", TABPFN_DEVICE, "\n\n",
  "NESTED-CV POOLED OOF\n",
  "LightGBM AUROC: ",
  sprintf("%.3f", lightgbm_auc),
  "\n",
  "TabPFN AUROC: ",
  sprintf("%.3f", tabpfn_auc),
  "\n",
  "Delta AUROC (TabPFN - LightGBM): ",
  sprintf("%+.3f", tabpfn_auc - lightgbm_auc),
  "\n",
  "Paired DeLong p: ",
  format(as.numeric(paired_delong$p.value), digits = 4),
  "\n\n",
  "STRICT LOCO\n",
  "LightGBM mean / worst AUROC: ",
  sprintf("%.3f", lgbm_loco$Mean_LOCO_AUROC[[1L]]),
  " / ",
  sprintf("%.3f", lgbm_loco$Worst_cohort_AUROC[[1L]]),
  "\n",
  "TabPFN mean / worst AUROC: ",
  sprintf("%.3f", tabpfn_loco$Mean_LOCO_AUROC[[1L]]),
  " / ",
  sprintf("%.3f", tabpfn_loco$Worst_cohort_AUROC[[1L]]),
  "\n\n",
  "DIAGNOSTIC FINAL ICE\n",
  "LightGBM AUROC: ",
  sprintf("%.3f", ice_lgbm_auc),
  "\n",
  "TabPFN AUROC: ",
  sprintf("%.3f", ice_tabpfn_auc),
  "\n",
  "Delta AUROC (TabPFN - LightGBM): ",
  sprintf("%+.3f", ice_tabpfn_auc - ice_lgbm_auc),
  "\n",
  "Paired DeLong p: ",
  format(as.numeric(ice_paired_delong$p.value), digits = 4),
  "\n\n",
  "Saved to:\n",
  FOUNDATION_DIR,
  "\n============================================================\n",
  sep = ""
)
