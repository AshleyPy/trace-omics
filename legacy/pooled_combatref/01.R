# ============================================================================
# CROSS-DISEASE CANDIDATE ADDITION
# VST + ComBat-Ref | Boruta top 50 | LightGBM
#
# PURPOSE
# -------
# Re-run the existing Model A / Model B development contrast using the
# study-derived candidate architecture, while preserving the existing 4,450-gene
# universe and the original development populations:
#
#   Model A: RA + UC + CD + T1D       (TB-naive; n = 1,413)
#   Model B: TB + RA + UC + CD + T1D  (ALL5;    n = 1,736)
#
# The script performs:
#   1. input and design audits;
#   2. leakage-free 5-fold x 2-repeat development CV;
#   3. fold-internal VST fitted on training data only;
#   4. fold-internal ComBat-Ref fitted on training cohorts only;
#   5. fold-internal Boruta top-50 feature selection;
#   6. fixed LightGBM fitting and OOF prediction;
#   7. development-only threshold selection (Youden) from averaged OOF predictions;
#   8. final refit on 100% development data;
#   9. independent shared final-ice evaluation WHEN all target cohorts are
#      represented in the relevant development population;
#  10. overall and disease-specific holdout metrics with uncertainty;
#  11. publication-ready audits and summary files.
#
# CRITICAL DESIGN RULE
# --------------------
# Frozen ComBat-Ref cannot transform a cohort/batch that was absent from the
# training data because neuroCombatFromTraining requires batch-specific
# estimates learned during training. Therefore:
#
#   - Model B can be evaluated on the shared final-ice set only if every final-ice
#     cohort is represented in Model B development.
#   - Model A will normally contain unseen TB cohorts in final ice. In STRICT_MODE
#     (default), the script WILL NOT estimate ComBat parameters from target data.
#     It saves an explicit blocker audit and skips that invalid final evaluation.
#
# This prevents target-data leakage / transductive harmonisation from being
# silently presented as strict cross-disease generalisability.
# ============================================================================

rm(list = ls(all.names = TRUE))
gc()

# ============================================================================
# 1. PACKAGES
# ============================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "SummarizedExperiment",
  "S4Vectors",
  "neuroCombat",
  "mlr3",
  "mlr3filters",
  "mlr3learners",
  "mlr3extralearners",
  "lightgbm",
  "pROC",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running:\n",
    paste(missing_packages, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(neuroCombat)
  library(mlr3)
  library(mlr3filters)
  library(mlr3learners)
  library(mlr3extralearners)
  library(lightgbm)
  library(pROC)
  library(ggplot2)
})

# ============================================================================
# 2. LOCKED SETTINGS
# ============================================================================

SEED <- 7832L
TOP_N <- 50L
CV_FOLDS <- 5L
CV_REPEATS <- 2L
BOOTSTRAPS <- 2000L
POSITIVE_CLASS <- "case"
NEGATIVE_CLASS <- "control"
STRICT_MODE <- TRUE
RESUME <- TRUE

# Exact fixed LightGBM specification used in the recent 7-model benchmark that
# generated the candidate architecture.
LGB_NUM_ITERATIONS <- 250L
LGB_MIN_DATA_IN_LEAF <- 5L
LGB_BAGGING_FRACTION <- 0.80
LGB_BAGGING_FREQ <- 1L

PIPELINE_LABEL <- "VST + ComBat-Ref | Boruta | None | LightGBM"

model_ids <- c(
  "Model_A_TB_naive",
  "Model_B_ALL5"
)

model_labels <- c(
  Model_A_TB_naive = "Model A: TB-naive",
  Model_B_ALL5 = "Model B: ALL5"
)

disease_order <- c("TB", "RA", "UC", "CD", "T1D")

set.seed(SEED)

# ============================================================================
# 3. PATHS
# ============================================================================

desktop_dir <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
root_dir <- file.path(desktop_dir, "Cross_disease_Two_Pooled_Models")

model_a_train_dir <- file.path(root_dir, "02_Model_A_TB_naive", "train")
model_b_train_dir <- file.path(root_dir, "03_Model_B_ALL5", "train")
shared_ice_dir <- file.path(root_dir, "04_SHARED_FINAL_ICE_TEST")

# NEW ADDITION. Nothing in the existing experiment is overwritten.
candidate_root <- file.path(
  root_dir,
  "12_CANDIDATE_FULL_VST_ComBatRef_Boruta_LightGBM"
)

fold_map_dir <- file.path(candidate_root, "01_fold_maps")
cv_preprocess_dir <- file.path(candidate_root, "02_CV_VST_ComBatRef")
cv_feature_dir <- file.path(candidate_root, "03_CV_Boruta")
cv_prediction_dir <- file.path(candidate_root, "04_CV_predictions")
cv_result_dir <- file.path(candidate_root, "05_CV_results")
threshold_dir <- file.path(candidate_root, "06_frozen_thresholds")
final_preprocess_dir <- file.path(candidate_root, "07_final_preprocessing")
final_feature_dir <- file.path(candidate_root, "08_final_Boruta")
final_model_dir <- file.path(candidate_root, "09_final_models")
final_prediction_dir <- file.path(candidate_root, "10_final_predictions")
final_result_dir <- file.path(candidate_root, "11_final_results")
audit_dir <- file.path(candidate_root, "12_audits")
figure_dir <- file.path(candidate_root, "13_figures")

for (d in c(
  candidate_root, fold_map_dir, cv_preprocess_dir, cv_feature_dir,
  cv_prediction_dir, cv_result_dir, threshold_dir, final_preprocess_dir,
  final_feature_dir, final_model_dir, final_prediction_dir,
  final_result_dir, audit_dir, figure_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# 4. INPUT FILES
# ============================================================================

development_paths <- list(
  Model_A_TB_naive = list(
    counts = file.path(
      model_a_train_dir,
      "Model_A_TB_naive_TRAIN_raw_counts_4450genes.rds"
    ),
    metadata = file.path(
      model_a_train_dir,
      "Model_A_TB_naive_TRAIN_metadata.rds"
    )
  ),
  Model_B_ALL5 = list(
    counts = file.path(
      model_b_train_dir,
      "Model_B_ALL5_TRAIN_raw_counts_4450genes.rds"
    ),
    metadata = file.path(
      model_b_train_dir,
      "Model_B_ALL5_TRAIN_metadata.rds"
    )
  )
)

resolve_existing_file <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(label, " was not found. Expected one of:\n", paste(candidates, collapse = "\n"))
  }
  existing[[1L]]
}

final_ice_counts_file <- resolve_existing_file(
  c(
    file.path(shared_ice_dir, "SHARED_ALL5_FINAL_ICE_TEST_raw_counts_4450genes.rds"),
    file.path(shared_ice_dir, "SHARED_ALL5_FINAL_ICE_TEST_raw_counts.rds")
  ),
  "Shared final-ice raw-count matrix"
)

final_ice_metadata_file <- resolve_existing_file(
  c(
    file.path(shared_ice_dir, "SHARED_ALL5_FINAL_ICE_TEST_metadata.rds"),
    file.path(shared_ice_dir, "SHARED_FINAL_ICE_TEST_metadata.rds")
  ),
  "Shared final-ice metadata"
)

# ============================================================================
# 5. GENERAL HELPERS
# ============================================================================

seed_for <- function(...) {
  vals <- as.integer(unlist(list(...)))
  as.integer((SEED + sum(vals * seq_along(vals) * 1009L)) %% .Machine$integer.max)
}

assert_finite_matrix <- function(x, label) {
  if (anyNA(x) || !all(is.finite(x))) {
    stop(label, " contains NA/Inf values.")
  }
  invisible(TRUE)
}

prepare_metadata <- function(metadata, expression_matrix, label) {
  metadata <- data.table::as.data.table(data.table::copy(metadata))
  expression_matrix <- as.matrix(expression_matrix)

  sample_candidates <- c("sample_id", "geo_accession", "sample")
  sample_column <- sample_candidates[sample_candidates %in% names(metadata)]
  if (length(sample_column) == 0L) {
    stop(label, ": no sample identifier was found.")
  }
  sample_column <- sample_column[[1L]]

  metadata[, sample_id := trimws(as.character(get(sample_column)))]

  if (!"outcome" %in% names(metadata)) stop(label, ": outcome column absent.")
  metadata[, outcome := tolower(trimws(as.character(outcome)))]
  if (!all(metadata$outcome %in% c(NEGATIVE_CLASS, POSITIVE_CLASS))) {
    stop(label, ": unexpected outcome labels.")
  }
  metadata[, outcome := factor(outcome, levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))]

  if (!"disease" %in% names(metadata)) stop(label, ": disease column absent.")
  metadata[, disease := trimws(as.character(disease))]

  if ("cohort_group" %in% names(metadata)) {
    metadata[, cohort_group := trimws(as.character(cohort_group))]
  } else if ("cohort_orig" %in% names(metadata)) {
    metadata[, cohort_group := paste(disease, trimws(as.character(cohort_orig)), sep = "__")]
  } else if ("cohort" %in% names(metadata)) {
    metadata[, cohort_group := paste(disease, trimws(as.character(cohort)), sep = "__")]
  } else {
    stop(label, ": cohort information absent.")
  }

  if (anyDuplicated(metadata$sample_id)) stop(label, ": duplicated metadata sample IDs.")
  if (is.null(colnames(expression_matrix))) stop(label, ": count matrix has no sample IDs.")

  idx <- match(colnames(expression_matrix), metadata$sample_id)
  if (anyNA(idx)) stop(label, ": count-matrix samples could not all be matched to metadata.")
  metadata <- metadata[idx]

  stopifnot(identical(colnames(expression_matrix), metadata$sample_id))
  metadata
}

validate_raw_counts <- function(counts, label) {
  counts <- as.matrix(counts)
  if (anyNA(counts) || !all(is.finite(counts)) || any(counts < 0)) {
    stop(label, ": invalid raw counts.")
  }
  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop(label, ": non-integer values found in raw counts.")
  }
  storage.mode(counts) <- "integer"
  counts
}

# ============================================================================
# 6. LOAD AND AUDIT DEVELOPMENT DATA
# ============================================================================

development_data <- list()

for (model_id in model_ids) {
  p <- development_paths[[model_id]]
  if (!file.exists(p$counts)) stop("Missing counts: ", p$counts)
  if (!file.exists(p$metadata)) stop("Missing metadata: ", p$metadata)

  counts <- validate_raw_counts(readRDS(p$counts), model_id)
  metadata <- prepare_metadata(readRDS(p$metadata), counts, model_id)

  development_data[[model_id]] <- list(
    counts = counts,
    metadata = metadata
  )
}

stopifnot(
  ncol(development_data[["Model_A_TB_naive"]]$counts) == 1413L,
  ncol(development_data[["Model_B_ALL5"]]$counts) == 1736L,
  nrow(development_data[["Model_A_TB_naive"]]$counts) == 4450L,
  nrow(development_data[["Model_B_ALL5"]]$counts) == 4450L
)

if ("TB" %in% development_data[["Model_A_TB_naive"]]$metadata$disease) {
  stop("TB has entered Model A development data.")
}

if (!setequal(unique(development_data[["Model_B_ALL5"]]$metadata$disease), disease_order)) {
  stop("Model B does not contain all five development diseases.")
}

reference_gene_order <- rownames(development_data[["Model_A_TB_naive"]]$counts)
if (!setequal(reference_gene_order, rownames(development_data[["Model_B_ALL5"]]$counts))) {
  stop("Model A and Model B do not share the same 4,450-gene universe.")
}

development_data[["Model_B_ALL5"]]$counts <-
  development_data[["Model_B_ALL5"]]$counts[reference_gene_order, , drop = FALSE]

# ============================================================================
# 7. LOAD SHARED FINAL-ICE DATA
# ============================================================================

final_ice_counts <- validate_raw_counts(readRDS(final_ice_counts_file), "final ice")
final_ice_metadata <- prepare_metadata(
  readRDS(final_ice_metadata_file),
  final_ice_counts,
  "final ice"
)

missing_genes <- setdiff(reference_gene_order, rownames(final_ice_counts))
if (length(missing_genes) > 0L) {
  stop("Final ice is missing genes from the fixed 4,450-gene universe: ",
       paste(head(missing_genes, 20L), collapse = ", "))
}
final_ice_counts <- final_ice_counts[reference_gene_order, , drop = FALSE]

# Design audit
input_audit <- rbindlist(lapply(model_ids, function(model_id) {
  md <- development_data[[model_id]]$metadata
  data.table(
    Model = model_labels[[model_id]],
    Development_n = nrow(md),
    Genes = nrow(development_data[[model_id]]$counts),
    Controls = sum(md$outcome == NEGATIVE_CLASS),
    Cases = sum(md$outcome == POSITIVE_CLASS),
    IR = sum(md$outcome == NEGATIVE_CLASS) / sum(md$outcome == POSITIVE_CLASS),
    Diseases = paste(sort(unique(md$disease)), collapse = " + "),
    Cohorts = uniqueN(md$cohort_group)
  )
}))

fwrite(input_audit, file.path(audit_dir, "01_development_input_audit.csv"))

final_ice_audit <- data.table(
  Evaluation_n = nrow(final_ice_metadata),
  Genes = nrow(final_ice_counts),
  Controls = sum(final_ice_metadata$outcome == NEGATIVE_CLASS),
  Cases = sum(final_ice_metadata$outcome == POSITIVE_CLASS),
  IR = sum(final_ice_metadata$outcome == NEGATIVE_CLASS) /
    sum(final_ice_metadata$outcome == POSITIVE_CLASS),
  Diseases = paste(sort(unique(final_ice_metadata$disease)), collapse = " + "),
  Cohorts = uniqueN(final_ice_metadata$cohort_group)
)

fwrite(final_ice_audit, file.path(audit_dir, "02_final_ice_input_audit.csv"))

# ============================================================================
# 8. REPEATED CV FOLD MAPS: COHORT x OUTCOME STRATIFIED
# ============================================================================

make_stratified_fold_map <- function(metadata, folds, repeats, seed) {
  rows <- list()
  k <- 0L

  for (repeat_index in seq_len(repeats)) {
    strata <- interaction(
      metadata$cohort_group,
      metadata$outcome,
      drop = TRUE,
      lex.order = TRUE
    )

    fold_assignment <- integer(nrow(metadata))
    stratum_levels <- levels(strata)

    for (s in seq_along(stratum_levels)) {
      stratum_name <- stratum_levels[[s]]
      idx <- which(strata == stratum_name)
      set.seed(seed + repeat_index * 10007L + s * 101L)
      idx <- sample(idx, length(idx), replace = FALSE)
      fold_assignment[idx] <- rep(seq_len(folds), length.out = length(idx))
    }

    for (i in seq_len(nrow(metadata))) {
      k <- k + 1L
      rows[[k]] <- data.table(
        sample_id = metadata$sample_id[[i]],
        disease = metadata$disease[[i]],
        cohort_group = metadata$cohort_group[[i]],
        outcome = as.character(metadata$outcome[[i]]),
        Repeat = repeat_index,
        Fold = fold_assignment[[i]]
      )
    }
  }

  out <- rbindlist(rows)

  # Every validation cohort must remain represented in its training partition.
  for (r in seq_len(repeats)) {
    for (f in seq_len(folds)) {
      val_ids <- out[Repeat == r & Fold == f, sample_id]
      train_md <- metadata[!sample_id %in% val_ids]
      val_md <- metadata[sample_id %in% val_ids]

      if (uniqueN(train_md$outcome) < 2L) {
        stop("CV training fold has only one outcome class.")
      }

      unseen <- setdiff(unique(val_md$cohort_group), unique(train_md$cohort_group))
      if (length(unseen) > 0L) {
        stop(
          "At least one validation cohort is absent from its CV training fold: ",
          paste(unseen, collapse = ", "),
          ". Reduce CV_FOLDS if necessary."
        )
      }
    }
  }

  out
}

fold_maps <- list()
for (model_id in model_ids) {
  fold_file <- file.path(fold_map_dir, paste0(model_id, "_5fold_2repeat.csv"))
  if (RESUME && file.exists(fold_file)) {
    fold_maps[[model_id]] <- fread(fold_file)
  } else {
    fold_maps[[model_id]] <- make_stratified_fold_map(
      development_data[[model_id]]$metadata,
      CV_FOLDS,
      CV_REPEATS,
      SEED
    )
    fwrite(fold_maps[[model_id]], fold_file)
  }
}

# ============================================================================
# 9. FROZEN VST
# ============================================================================

calculate_poscounts_geo_means <- function(count_matrix) {
  positive_n <- rowSums(count_matrix > 0)
  if (any(positive_n == 0L)) {
    stop("At least one retained gene is zero in every training observation.")
  }

  log_counts <- log(count_matrix)
  log_counts[!is.finite(log_counts)] <- 0
  exp(rowSums(log_counts) / positive_n)
}

fit_apply_frozen_vst <- function(train_counts, new_counts) {
  common_genes <- intersect(rownames(train_counts), rownames(new_counts))
  train_counts <- train_counts[common_genes, , drop = FALSE]
  new_counts <- new_counts[common_genes, , drop = FALSE]

  geo_means <- calculate_poscounts_geo_means(train_counts)

  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData = train_counts,
    colData = S4Vectors::DataFrame(row.names = colnames(train_counts)),
    design = ~ 1
  )

  dds_train <- DESeq2::estimateSizeFactors(dds_train, geoMeans = geo_means)

  dispersion_error <- NULL
  for (fit_type in c("parametric", "local", "mean")) {
    tmp <- tryCatch(
      suppressMessages(DESeq2::estimateDispersions(dds_train, fitType = fit_type)),
      error = function(e) {
        dispersion_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(tmp)) {
      dds_train <- tmp
      break
    }
  }

  if (is.null(DESeq2::dispersionFunction(dds_train))) {
    stop("DESeq2 dispersion estimation failed: ", dispersion_error)
  }

  train_dispersion_function <- DESeq2::dispersionFunction(dds_train)

  vsd_train <- suppressMessages(
    DESeq2::varianceStabilizingTransformation(dds_train, blind = FALSE)
  )
  expr_train <- SummarizedExperiment::assay(vsd_train)

  dds_new <- DESeq2::DESeqDataSetFromMatrix(
    countData = new_counts,
    colData = S4Vectors::DataFrame(row.names = colnames(new_counts)),
    design = ~ 1
  )
  dds_new <- DESeq2::estimateSizeFactors(dds_new, geoMeans = geo_means)
  DESeq2::dispersionFunction(dds_new) <- train_dispersion_function

  vsd_new <- suppressMessages(
    DESeq2::varianceStabilizingTransformation(dds_new, blind = FALSE)
  )
  expr_new <- SummarizedExperiment::assay(vsd_new)

  dimnames(expr_train) <- dimnames(train_counts)
  dimnames(expr_new) <- dimnames(new_counts)

  assert_finite_matrix(expr_train, "Frozen VST training matrix")
  assert_finite_matrix(expr_new, "Frozen VST new matrix")

  list(
    train = expr_train,
    new = expr_new,
    state = list(
      genes = rownames(expr_train),
      geo_means = geo_means,
      dispersion_function = train_dispersion_function
    )
  )
}

# ============================================================================
# 10. FROZEN COMBAT-REF
# ============================================================================

within_batch_gene_filter <- function(expr_train, batch_train) {
  residual_ss <- numeric(nrow(expr_train))
  for (b in unique(batch_train)) {
    cols <- which(batch_train == b)
    current <- expr_train[, cols, drop = FALSE]
    centred <- sweep(current, 1L, rowMeans(current), FUN = "-")
    residual_ss <- residual_ss + rowSums(centred^2)
  }
  is.finite(residual_ss) & residual_ss > sqrt(.Machine$double.eps)
}

apply_neurocombat_from_training <- function(new_dat, new_batch, estimates) {
  tryCatch(
    neuroCombat::neuroCombatFromTraining(
      dat = new_dat,
      batch = new_batch,
      estimates = estimates,
      mod = NULL,
      verbose = FALSE
    ),
    error = function(e1) {
      neuroCombat::neuroCombatFromTraining(
        dat = new_dat,
        batch = new_batch,
        estimates = estimates,
        verbose = FALSE
      )
    }
  )
}

fit_apply_frozen_combat_ref <- function(expr_train, expr_new, batch_train, batch_new) {
  batch_train <- as.character(batch_train)
  batch_new <- as.character(batch_new)

  unseen_batches <- setdiff(unique(batch_new), unique(batch_train))
  if (length(unseen_batches) > 0L) {
    stop(
      "STRICT_COMBAT_UNSEEN_BATCH: ",
      paste(unseen_batches, collapse = ", ")
    )
  }

  keep <- within_batch_gene_filter(expr_train, batch_train)
  if (!any(keep)) stop("No genes had sufficient within-cohort training variation for ComBat-Ref.")

  train_input <- expr_train[keep, , drop = FALSE]
  new_input <- expr_new[keep, , drop = FALSE]

  ref_level <- names(sort(table(batch_train), decreasing = TRUE))[1L]

  fit <- neuroCombat::neuroCombat(
    dat = train_input,
    batch = batch_train,
    mod = NULL,
    eb = TRUE,
    parametric = TRUE,
    mean.only = FALSE,
    ref.batch = ref_level,
    verbose = FALSE
  )

  if (is.null(fit$dat.combat) || is.null(fit$estimates)) {
    stop("neuroCombat did not return adjusted data and reusable estimates.")
  }

  new_fit <- apply_neurocombat_from_training(
    new_dat = new_input,
    new_batch = batch_new,
    estimates = fit$estimates
  )

  if (is.null(new_fit$dat.combat)) {
    stop("neuroCombatFromTraining did not return adjusted new data.")
  }

  adjusted_train <- expr_train
  adjusted_new <- expr_new
  adjusted_train[keep, ] <- as.matrix(fit$dat.combat)
  adjusted_new[keep, ] <- as.matrix(new_fit$dat.combat)

  dimnames(adjusted_train) <- dimnames(expr_train)
  dimnames(adjusted_new) <- dimnames(expr_new)

  assert_finite_matrix(adjusted_train, "ComBat-Ref training matrix")
  assert_finite_matrix(adjusted_new, "ComBat-Ref new matrix")

  list(
    train = adjusted_train,
    new = adjusted_new,
    state = list(
      reference_cohort = ref_level,
      estimates = fit$estimates,
      genes_estimated = rownames(expr_train)[keep],
      genes_unchanged = rownames(expr_train)[!keep]
    )
  )
}

# ============================================================================
# 11. BORUTA TOP-50 FEATURE SELECTION
# ============================================================================

make_feature_task <- function(expr_train, metadata, task_id) {
  original_genes <- rownames(expr_train)
  safe_genes <- make.names(original_genes, unique = TRUE)

  x <- data.table::as.data.table(t(expr_train))
  data.table::setnames(x, safe_genes)
  x[, outcome := factor(as.character(metadata$outcome),
                        levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))]

  task <- mlr3::TaskClassif$new(
    id = task_id,
    backend = x,
    target = "outcome",
    positive = POSITIVE_CLASS
  )

  map <- data.table(
    safe_gene = safe_genes,
    original_gene = original_genes
  )

  list(task = task, map = map)
}

normalise_filter_output <- function(result, filter_object) {
  if (inherits(result, "data.frame") || data.table::is.data.table(result)) {
    dt <- data.table::as.data.table(result)
    feature_col <- c("feature", "feature_name", "name")
    feature_col <- feature_col[feature_col %in% names(dt)]
    score_col <- c("score", "importance", "value")
    score_col <- score_col[score_col %in% names(dt)]
    if (length(feature_col) > 0L && length(score_col) > 0L) {
      return(dt[, .(
        safe_gene = as.character(get(feature_col[[1L]])),
        score = as.numeric(get(score_col[[1L]]))
      )])
    }
  }

  scores <- filter_object$scores
  if (is.null(scores) || is.null(names(scores))) {
    stop("Boruta returned an unsupported score structure.")
  }

  data.table(
    safe_gene = names(scores),
    score = as.numeric(scores)
  )
}

select_boruta_top50 <- function(expr_train, metadata, seed, task_id) {
  task_info <- make_feature_task(expr_train, metadata, task_id)
  set.seed(seed)
  filter_object <- mlr3filters::flt("boruta")
  result <- filter_object$calculate(task_info$task)
  score_dt <- normalise_filter_output(result, filter_object)

  score_dt <- score_dt[
    !is.na(safe_gene) & nzchar(safe_gene) & is.finite(score)
  ]

  score_dt[, original_gene := task_info$map$original_gene[
    match(safe_gene, task_info$map$safe_gene)
  ]]
  score_dt <- score_dt[!is.na(original_gene)]
  setorder(score_dt, -score, safe_gene)
  score_dt <- unique(score_dt, by = "original_gene")

  if (nrow(score_dt) < TOP_N) {
    stop("Boruta returned fewer than ", TOP_N, " finite ranked genes.")
  }

  score_dt[, Rank := seq_len(.N)]
  selected <- score_dt[seq_len(TOP_N), original_gene]

  list(genes = selected, scores = score_dt)
}

# ============================================================================
# 12. LIGHTGBM
# ============================================================================

make_fixed_lightgbm <- function() {
  if (!"classif.lightgbm" %in% mlr3::mlr_learners$keys()) {
    stop("classif.lightgbm is not registered in mlr3.")
  }

  mlr3::lrn(
    "classif.lightgbm",
    id = "LightGBM",
    predict_type = "prob",
    num_iterations = LGB_NUM_ITERATIONS,
    min_data_in_leaf = LGB_MIN_DATA_IN_LEAF,
    bagging_fraction = LGB_BAGGING_FRACTION,
    bagging_freq = LGB_BAGGING_FREQ,
    num_threads = 1L,
    device_type = "cpu",
    deterministic = TRUE
  )
}

make_model_task <- function(expr, metadata, task_id) {
  original_genes <- rownames(expr)
  safe_genes <- make.names(original_genes, unique = TRUE)
  x <- data.table::as.data.table(t(expr))
  data.table::setnames(x, safe_genes)
  x[, outcome := factor(as.character(metadata$outcome),
                        levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))]

  mlr3::TaskClassif$new(
    id = task_id,
    backend = x,
    target = "outcome",
    positive = POSITIVE_CLASS
  )
}

# ============================================================================
# 13. METRICS
# ============================================================================

safe_roc <- function(truth, probability) {
  truth <- factor(as.character(truth), levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
  if (length(unique(truth)) < 2L) return(NULL)
  pROC::roc(
    response = truth,
    predictor = probability,
    levels = c(NEGATIVE_CLASS, POSITIVE_CLASS),
    direction = "<",
    quiet = TRUE
  )
}

safe_auroc <- function(truth, probability) {
  r <- safe_roc(truth, probability)
  if (is.null(r)) return(NA_real_)
  as.numeric(pROC::auc(r))
}

# Average precision / step-wise PR area.
safe_auprc <- function(truth, probability) {
  y <- as.integer(as.character(truth) == POSITIVE_CLASS)
  ok <- is.finite(probability) & !is.na(y)
  y <- y[ok]
  p <- probability[ok]
  if (sum(y == 1L) == 0L) return(NA_real_)

  ord <- order(p, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / (tp + fp)
  recall <- tp / sum(y == 1L)
  recall_prev <- c(0, head(recall, -1L))
  sum((recall - recall_prev) * precision)
}

metric_vector <- function(truth, probability, threshold = 0.50) {
  truth_chr <- as.character(truth)
  pred <- ifelse(probability >= threshold, POSITIVE_CLASS, NEGATIVE_CLASS)

  tp <- sum(pred == POSITIVE_CLASS & truth_chr == POSITIVE_CLASS)
  tn <- sum(pred == NEGATIVE_CLASS & truth_chr == NEGATIVE_CLASS)
  fp <- sum(pred == POSITIVE_CLASS & truth_chr == NEGATIVE_CLASS)
  fn <- sum(pred == NEGATIVE_CLASS & truth_chr == POSITIVE_CLASS)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

  c(
    AUROC = safe_auroc(truth, probability),
    AUPRC = safe_auprc(truth, probability),
    BACC = mean(c(sensitivity, specificity), na.rm = TRUE),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Brier = mean((probability - as.integer(truth_chr == POSITIVE_CLASS))^2)
  )
}

choose_youden_threshold <- function(truth, probability) {
  roc_obj <- safe_roc(truth, probability)
  if (is.null(roc_obj)) stop("Cannot select threshold from one-class OOF predictions.")

  coords <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )

  threshold <- as.numeric(coords$threshold[[1L]])
  if (!is.finite(threshold)) stop("Youden threshold was non-finite.")
  threshold
}

bootstrap_metric_ci <- function(truth, probability, threshold, B = BOOTSTRAPS, seed = SEED) {
  truth_chr <- as.character(truth)
  case_idx <- which(truth_chr == POSITIVE_CLASS)
  control_idx <- which(truth_chr == NEGATIVE_CLASS)

  point <- metric_vector(truth, probability, threshold)
  boot <- matrix(NA_real_, nrow = B, ncol = length(point))
  colnames(boot) <- names(point)

  set.seed(seed)
  for (b in seq_len(B)) {
    idx <- c(
      sample(case_idx, length(case_idx), replace = TRUE),
      sample(control_idx, length(control_idx), replace = TRUE)
    )
    boot[b, ] <- metric_vector(truth[idx], probability[idx], threshold)
  }

  out <- rbindlist(lapply(names(point), function(metric_name) {
    vals <- boot[, metric_name]
    data.table(
      Metric = metric_name,
      Estimate = unname(point[[metric_name]]),
      Lower95 = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
      Upper95 = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE))
    )
  }))

  # Replace AUROC bootstrap CI with DeLong CI when possible.
  roc_obj <- safe_roc(truth, probability)
  if (!is.null(roc_obj)) {
    ci_auc <- as.numeric(pROC::ci.auc(roc_obj, method = "delong", conf.level = 0.95))
    out[Metric == "AUROC", `:=`(Lower95 = ci_auc[[1L]], Upper95 = ci_auc[[3L]])]
  }

  out
}

# ============================================================================
# 14. ONE LEAKAGE-FREE CV FOLD
# ============================================================================

run_one_cv_fold <- function(model_id, repeat_index, fold_index) {
  cache_file <- file.path(
    cv_result_dir,
    paste0(model_id, "__R", repeat_index, "__F", fold_index, "__candidate.rds")
  )
  if (RESUME && file.exists(cache_file)) return(readRDS(cache_file))

  dat <- development_data[[model_id]]
  fold_map <- fold_maps[[model_id]]

  val_ids <- fold_map[Repeat == repeat_index & Fold == fold_index, sample_id]
  train_ids <- setdiff(dat$metadata$sample_id, val_ids)

  train_idx <- match(train_ids, dat$metadata$sample_id)
  val_idx <- match(val_ids, dat$metadata$sample_id)

  train_counts <- dat$counts[, train_idx, drop = FALSE]
  val_counts <- dat$counts[, val_idx, drop = FALSE]
  train_md <- dat$metadata[train_idx]
  val_md <- dat$metadata[val_idx]

  stopifnot(
    identical(colnames(train_counts), train_md$sample_id),
    identical(colnames(val_counts), val_md$sample_id)
  )

  vst <- fit_apply_frozen_vst(train_counts, val_counts)

  combat <- fit_apply_frozen_combat_ref(
    vst$train,
    vst$new,
    train_md$cohort_group,
    val_md$cohort_group
  )

  fs <- select_boruta_top50(
    combat$train,
    train_md,
    seed = seed_for(repeat_index, fold_index, match(model_id, model_ids)),
    task_id = paste0("Boruta_", model_id, "_R", repeat_index, "_F", fold_index)
  )

  train_expr <- combat$train[fs$genes, , drop = FALSE]
  val_expr <- combat$new[fs$genes, , drop = FALSE]

  train_task <- make_model_task(
    train_expr,
    train_md,
    paste0("train_", model_id, "_R", repeat_index, "_F", fold_index)
  )
  val_task <- make_model_task(
    val_expr,
    val_md,
    paste0("val_", model_id, "_R", repeat_index, "_F", fold_index)
  )

  learner <- make_fixed_lightgbm()
  set.seed(seed_for(900L, repeat_index, fold_index, match(model_id, model_ids)))
  learner$train(train_task)
  pred <- learner$predict(val_task)
  prob <- pred$prob[, POSITIVE_CLASS]

  metrics_05 <- metric_vector(val_md$outcome, prob, threshold = 0.50)

  pred_dt <- data.table(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Pipeline = PIPELINE_LABEL,
    Repeat = repeat_index,
    Fold = fold_index,
    sample_id = val_md$sample_id,
    disease = val_md$disease,
    cohort_group = val_md$cohort_group,
    truth = as.character(val_md$outcome),
    probability_case = prob
  )

  feature_dt <- fs$scores[, .(
    Model_ID = model_id,
    Repeat = repeat_index,
    Fold = fold_index,
    Rank,
    Gene = original_gene,
    Score = score,
    Selected = original_gene %in% fs$genes
  )]

  result <- list(
    prediction = pred_dt,
    metrics_05 = metrics_05,
    selected_genes = fs$genes,
    feature_scores = feature_dt,
    reference_cohort = combat$state$reference_cohort,
    train_n = nrow(train_md),
    validation_n = nrow(val_md)
  )

  fwrite(
    pred_dt,
    file.path(cv_prediction_dir,
              paste0(model_id, "__R", repeat_index, "__F", fold_index, "__predictions.csv"))
  )
  fwrite(
    feature_dt,
    file.path(cv_feature_dir,
              paste0(model_id, "__R", repeat_index, "__F", fold_index, "__Boruta_scores.csv"))
  )
  saveRDS(
    list(vst_state = vst$state, combat_state = combat$state),
    file.path(cv_preprocess_dir,
              paste0(model_id, "__R", repeat_index, "__F", fold_index, "__preprocessing_state.rds"))
  )
  saveRDS(result, cache_file)
  result
}

# ============================================================================
# 15. RUN DEVELOPMENT CV
# ============================================================================

cv_results <- list()
cv_fold_rows <- list()
cv_prediction_rows <- list()
cv_feature_rows <- list()
job <- 0L

for (model_id in model_ids) {
  for (r in seq_len(CV_REPEATS)) {
    for (f in seq_len(CV_FOLDS)) {
      job <- job + 1L
      cat("\nRunning candidate CV: ", model_id, " repeat ", r, " fold ", f, "\n", sep = "")
      res <- run_one_cv_fold(model_id, r, f)
      cv_results[[paste(model_id, r, f, sep = "__")]] <- res

      cv_prediction_rows[[job]] <- res$prediction
      cv_feature_rows[[job]] <- res$feature_scores[Selected == TRUE]

      cv_fold_rows[[job]] <- data.table(
        Model_ID = model_id,
        Model = model_labels[[model_id]],
        Repeat = r,
        Fold = f,
        AUROC = res$metrics_05[["AUROC"]],
        AUPRC = res$metrics_05[["AUPRC"]],
        BACC = res$metrics_05[["BACC"]],
        Sensitivity = res$metrics_05[["Sensitivity"]],
        Specificity = res$metrics_05[["Specificity"]],
        Brier = res$metrics_05[["Brier"]],
        Reference_cohort = res$reference_cohort,
        Train_n = res$train_n,
        Validation_n = res$validation_n
      )
    }
  }
}

cv_fold_metrics <- rbindlist(cv_fold_rows)
cv_predictions <- rbindlist(cv_prediction_rows)
cv_selected_features <- rbindlist(cv_feature_rows)

fwrite(cv_fold_metrics, file.path(cv_result_dir, "ALL_candidate_CV_fold_metrics.csv"))
fwrite(cv_predictions, file.path(cv_prediction_dir, "ALL_candidate_CV_predictions.csv"))
fwrite(cv_selected_features, file.path(cv_feature_dir, "ALL_candidate_CV_selected_genes.csv"))

cv_summary <- cv_fold_metrics[, .(
  Folds = .N,
  Mean_AUROC = mean(AUROC, na.rm = TRUE),
  SD_AUROC = sd(AUROC, na.rm = TRUE),
  Mean_AUPRC = mean(AUPRC, na.rm = TRUE),
  Mean_BACC_0.50 = mean(BACC, na.rm = TRUE),
  Mean_Sensitivity_0.50 = mean(Sensitivity, na.rm = TRUE),
  Mean_Specificity_0.50 = mean(Specificity, na.rm = TRUE),
  Mean_Brier = mean(Brier, na.rm = TRUE)
), by = .(Model_ID, Model)]

fwrite(cv_summary, file.path(cv_result_dir, "candidate_CV_summary.csv"))

# ============================================================================
# 16. FROZEN THRESHOLDS FROM DEVELOPMENT OOF ONLY
# ============================================================================

frozen_thresholds <- list()
threshold_rows <- list()

for (model_id in model_ids) {
  pred_model <- cv_predictions[Model_ID == model_id]

  # Each sample appears twice because CV_REPEATS = 2. Average its two OOF
  # probabilities before selecting the operating threshold.
  oof_avg <- pred_model[, .(
    truth = unique(truth)[1L],
    probability_case = mean(probability_case, na.rm = TRUE)
  ), by = .(sample_id)]

  threshold <- choose_youden_threshold(oof_avg$truth, oof_avg$probability_case)
  frozen_thresholds[[model_id]] <- threshold

  m <- metric_vector(oof_avg$truth, oof_avg$probability_case, threshold)
  threshold_rows[[model_id]] <- data.table(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Threshold_source = "Averaged repeated-CV OOF predictions only",
    Frozen_threshold = threshold,
    OOF_AUROC = m[["AUROC"]],
    OOF_AUPRC = m[["AUPRC"]],
    OOF_BACC = m[["BACC"]],
    OOF_Sensitivity = m[["Sensitivity"]],
    OOF_Specificity = m[["Specificity"]],
    OOF_Brier = m[["Brier"]]
  )

  saveRDS(
    list(threshold = threshold, oof_average = oof_avg),
    file.path(threshold_dir, paste0(model_id, "__frozen_threshold.rds"))
  )
}

threshold_table <- rbindlist(threshold_rows)
fwrite(threshold_table, file.path(threshold_dir, "candidate_frozen_thresholds.csv"))

# ============================================================================
# 17. FINAL FIT + STRICT FINAL-ICE EVALUATION
# ============================================================================

fit_final_candidate <- function(model_id) {
  dat <- development_data[[model_id]]
  threshold <- frozen_thresholds[[model_id]]

  train_batches <- unique(dat$metadata$cohort_group)
  target_batches <- unique(final_ice_metadata$cohort_group)
  unseen_batches <- setdiff(target_batches, train_batches)

  batch_audit <- data.table(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Training_cohorts = length(train_batches),
    Final_ice_cohorts = length(target_batches),
    Unseen_final_ice_cohorts = length(unseen_batches),
    Unseen_cohort_names = paste(unseen_batches, collapse = " | "),
    Strict_ComBatRef_possible = length(unseen_batches) == 0L
  )
  fwrite(batch_audit, file.path(audit_dir, paste0(model_id, "__ComBatRef_target_batch_audit.csv")))

  if (length(unseen_batches) > 0L && STRICT_MODE) {
    message(
      "\nSTRICT MODE: skipping final-ice candidate evaluation for ",
      model_labels[[model_id]],
      " because these target cohorts were absent from development:\n  ",
      paste(unseen_batches, collapse = "\n  "),
      "\nNo target-data ComBat estimation was performed."
    )

    return(list(
      status = "SKIPPED_STRICT_UNSEEN_BATCH",
      batch_audit = batch_audit,
      unseen_batches = unseen_batches
    ))
  }

  vst <- fit_apply_frozen_vst(dat$counts, final_ice_counts)

  combat <- fit_apply_frozen_combat_ref(
    vst$train,
    vst$new,
    dat$metadata$cohort_group,
    final_ice_metadata$cohort_group
  )

  fs <- select_boruta_top50(
    combat$train,
    dat$metadata,
    seed = seed_for(5000L, match(model_id, model_ids)),
    task_id = paste0("FINAL_Boruta_", model_id)
  )

  train_expr <- combat$train[fs$genes, , drop = FALSE]
  test_expr <- combat$new[fs$genes, , drop = FALSE]

  train_task <- make_model_task(train_expr, dat$metadata, paste0("FINAL_train_", model_id))
  test_task <- make_model_task(test_expr, final_ice_metadata, paste0("FINAL_ice_", model_id))

  learner <- make_fixed_lightgbm()
  set.seed(seed_for(6000L, match(model_id, model_ids)))
  learner$train(train_task)
  pred <- learner$predict(test_task)
  prob <- pred$prob[, POSITIVE_CLASS]

  pred_dt <- data.table(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Pipeline = PIPELINE_LABEL,
    sample_id = final_ice_metadata$sample_id,
    disease = final_ice_metadata$disease,
    cohort_group = final_ice_metadata$cohort_group,
    truth = as.character(final_ice_metadata$outcome),
    probability_case = prob,
    Frozen_threshold = threshold,
    predicted_class = ifelse(prob >= threshold, POSITIVE_CLASS, NEGATIVE_CLASS)
  )

  # Save final preprocessing state and selected genes.
  saveRDS(
    list(vst_state = vst$state, combat_state = combat$state),
    file.path(final_preprocess_dir, paste0(model_id, "__final_preprocessing_state.rds"))
  )
  fwrite(
    fs$scores,
    file.path(final_feature_dir, paste0(model_id, "__FINAL_Boruta_scores.csv"))
  )
  saveRDS(fs$genes, file.path(final_feature_dir, paste0(model_id, "__FINAL_Boruta_top50_genes.rds")))
  fwrite(data.table(Rank = seq_along(fs$genes), Gene = fs$genes),
         file.path(final_feature_dir, paste0(model_id, "__FINAL_Boruta_top50_genes.csv")))

  # Save mlr3 learner. Note: lightgbm versions can affect portability of serialized boosters.
  saveRDS(learner, file.path(final_model_dir, paste0(model_id, "__FINAL_LightGBM_learner.rds")))
  fwrite(pred_dt, file.path(final_prediction_dir, paste0(model_id, "__FINAL_ICE_predictions.csv")))

  overall_ci <- bootstrap_metric_ci(
    pred_dt$truth,
    pred_dt$probability_case,
    threshold,
    B = BOOTSTRAPS,
    seed = seed_for(7000L, match(model_id, model_ids))
  )
  overall_ci[, `:=`(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Scope = "Overall",
    Disease = "ALL"
  )]

  disease_ci <- rbindlist(lapply(disease_order, function(disease_name) {
    sub <- pred_dt[disease == disease_name]
    if (nrow(sub) == 0L || uniqueN(sub$truth) < 2L) return(NULL)

    ans <- bootstrap_metric_ci(
      sub$truth,
      sub$probability_case,
      threshold,
      B = BOOTSTRAPS,
      seed = seed_for(8000L, match(disease_name, disease_order), match(model_id, model_ids))
    )
    ans[, `:=`(
      Model_ID = model_id,
      Model = model_labels[[model_id]],
      Scope = "Disease-specific",
      Disease = disease_name
    )]
    ans
  }), fill = TRUE)

  all_metrics <- rbindlist(list(overall_ci, disease_ci), fill = TRUE)
  fwrite(all_metrics, file.path(final_result_dir, paste0(model_id, "__FINAL_ICE_metrics_95CI.csv")))

  list(
    status = "COMPLETE",
    predictions = pred_dt,
    metrics = all_metrics,
    selected_genes = fs$genes,
    reference_cohort = combat$state$reference_cohort,
    threshold = threshold,
    batch_audit = batch_audit
  )
}

final_results <- lapply(model_ids, fit_final_candidate)
names(final_results) <- model_ids

# ============================================================================
# 18. COMBINE FINAL RESULTS + EXPLICIT STRICTNESS AUDIT
# ============================================================================

final_status <- rbindlist(lapply(model_ids, function(model_id) {
  res <- final_results[[model_id]]
  data.table(
    Model_ID = model_id,
    Model = model_labels[[model_id]],
    Pipeline = PIPELINE_LABEL,
    Status = res$status,
    Strict_mode = STRICT_MODE,
    Reason = if (res$status == "COMPLETE") {
      "All final-ice cohorts represented in development"
    } else {
      paste0(
        "Frozen ComBat-Ref cannot transform unseen final-ice cohorts: ",
        paste(res$unseen_batches, collapse = " | ")
      )
    }
  )
}))

fwrite(final_status, file.path(audit_dir, "03_candidate_final_evaluation_status.csv"))

completed_ids <- final_status[Status == "COMPLETE", Model_ID]

if (length(completed_ids) > 0L) {
  all_final_predictions <- rbindlist(lapply(completed_ids, function(id) {
    final_results[[id]]$predictions
  }), fill = TRUE)

  all_final_metrics <- rbindlist(lapply(completed_ids, function(id) {
    final_results[[id]]$metrics
  }), fill = TRUE)

  fwrite(all_final_predictions, file.path(final_prediction_dir, "ALL_candidate_FINAL_ICE_predictions.csv"))
  fwrite(all_final_metrics, file.path(final_result_dir, "ALL_candidate_FINAL_ICE_metrics_95CI.csv"))
}

# ============================================================================
# 19. DEVELOPMENT CONTRAST FIGURE
# ============================================================================

p_dev <- ggplot(cv_fold_metrics, aes(x = Model, y = AUROC)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 1.8, alpha = 0.75) +
  geom_hline(yintercept = 0.50, linetype = 2) +
  labs(
    x = NULL,
    y = "Fold AUROC",
    title = "Candidate-pipeline development discrimination"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "Figure_candidate_ModelA_vs_ModelB_development_AUROC.pdf"),
  p_dev, width = 6.8, height = 4.8, units = "in", device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "Figure_candidate_ModelA_vs_ModelB_development_AUROC.png"),
  p_dev, width = 6.8, height = 4.8, units = "in", dpi = 600
)

# ============================================================================
# 20. FINAL MODEL-B DISEASE-SPECIFIC FIGURE, IF AVAILABLE
# ============================================================================

if ("Model_B_ALL5" %in% completed_ids) {
  b_metrics <- final_results[["Model_B_ALL5"]]$metrics[
    Scope == "Disease-specific" & Metric == "AUROC"
  ]
  b_metrics[, Disease := factor(Disease, levels = disease_order)]

  p_b <- ggplot(b_metrics, aes(x = Estimate, y = Disease)) +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.15) +
    geom_point(size = 2.4) +
    geom_vline(xintercept = 0.50, linetype = 2) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      x = "AUROC (95% CI)",
      y = NULL,
      title = "ALL5 candidate performance by disease"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(figure_dir, "Figure_candidate_ALL5_disease_specific_AUROC.pdf"),
    p_b, width = 7.0, height = 4.8, units = "in", device = cairo_pdf
  )

  ggsave(
    file.path(figure_dir, "Figure_candidate_ALL5_disease_specific_AUROC.png"),
    p_b, width = 7.0, height = 4.8, units = "in", dpi = 600
  )
}
