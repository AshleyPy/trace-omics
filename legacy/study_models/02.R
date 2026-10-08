# ======================================================================
# FIX CD OUTCOME STANDARDISATION
# ======================================================================

standardise_tb_outcome <- function(
  outcome
) {

  original <- trimws(
    as.character(
      outcome
    )
  )

  lower <- tolower(
    original
  )


  case_labels <- c(
    "case",
    "cases",
    "cd",
    "crohn",
    "crohn's disease",
    "crohns disease",
    "disease",
    "diseased",
    "positive"
  )


  control_labels <- c(
    "control",
    "controls",
    "healthy",
    "non-cd",
    "noncd",
    "non-ibd",
    "negative",
    "hc"
  )


  standardised <- rep(
    NA_character_,
    length(
      lower
    )
  )


  standardised[
    lower %in%
      case_labels
  ] <- POSITIVE_CLASS


  standardised[
    lower %in%
      control_labels
  ] <- NEGATIVE_CLASS


  invalid <- unique(
    original[
      is.na(
        standardised
      )
    ]
  )


  if (
    length(
      invalid
    ) >
      0L
  ) {

    stop(
      "Unexpected CD outcome label(s):\n",
      paste(
        invalid,
        collapse = "\n"
      )
    )
  }


  factor(
    standardised,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )
}

# ======================================================================
# LOAD + ALIGN FULL TB DATA
# ======================================================================

tb_counts <- as.matrix(readRDS(counts_path))
tb_metadata <- data.table::as.data.table(readRDS(metadata_path))
if (is.null(rownames(tb_counts)) || is.null(colnames(tb_counts))) stop("CD count matrix needs row and column names.")
storage.mode(tb_counts) <- "numeric"
if (anyNA(tb_counts) || !all(is.finite(tb_counts)) || any(tb_counts < 0)) stop("Invalid raw counts.")
if (any(abs(tb_counts - round(tb_counts)) > 1e-8)) stop("CD_joined_numeric_matrix.rds is not integer-valued raw counts.")
tb_counts <- round(tb_counts); storage.mode(tb_counts) <- "integer"
if (anyDuplicated(rownames(tb_counts)) || anyDuplicated(colnames(tb_counts))) stop("Duplicated count identifiers detected.")
if (!"geo_accession" %in% names(tb_metadata)) stop("CD metadata requires geo_accession.")
if (!"outcome" %in% names(tb_metadata)) stop("CD metadata requires outcome.")
if ("cohort_group" %in% names(tb_metadata)) {
  tb_metadata[, cohort_group := trimws(as.character(cohort_group))]
  tb_metadata[, cohort_group := sub("^CD__", "", cohort_group)]
} else if ("cohort" %in% names(tb_metadata)) {
  tb_metadata[, cohort_group := trimws(as.character(cohort))]
} else stop("CD metadata requires cohort_group or cohort.")
tb_metadata[, sample_id := trimws(as.character(geo_accession))]
tb_metadata[, outcome := standardise_tb_outcome(outcome)]
if (anyDuplicated(tb_metadata$sample_id)) stop("Duplicated CD metadata sample IDs.")
common_samples <- colnames(tb_counts)[colnames(tb_counts) %in% tb_metadata$sample_id]
tb_counts <- tb_counts[, common_samples, drop = FALSE]
tb_metadata <- tb_metadata[match(common_samples, sample_id)]
stopifnot(identical(colnames(tb_counts), tb_metadata$sample_id))
tb_metadata[, disease := "CD"]

if (nrow(tb_metadata) != 495L) {
  warning("Expected 495 aligned CD samples; observed ", nrow(tb_metadata), ". The script will use the observed aligned population.")
}

raw_audit <- data.table::data.table(
  Population = "Full CD", n = nrow(tb_metadata),
  Controls = sum(tb_metadata$outcome == NEGATIVE_CLASS),
  Cases = sum(tb_metadata$outcome == POSITIVE_CLASS),
  Cohorts = data.table::uniqueN(tb_metadata$cohort_group), Genes = nrow(tb_counts)
)
data.table::fwrite(raw_audit, file.path(audit_dir, "CD_full_raw_input_audit.csv"))
print(raw_audit)
print(table(tb_metadata$cohort_group, tb_metadata$outcome))

# ======================================================================
# EXACT APPROXIMATELY 80/20 COHORT x OUTCOME-STRATIFIED FINAL ICE SPLIT
# ======================================================================

allocate_integer_target <- function(counts_by_stratum, target, retain_at_least_one = TRUE) {
  counts_by_stratum <- as.integer(counts_by_stratum) |> stats::setNames(names(counts_by_stratum))
  eligible <- counts_by_stratum > 0L
  capacity <- if (retain_at_least_one) pmax(counts_by_stratum - as.integer(eligible), 0L) else counts_by_stratum
  if (sum(capacity) < target) stop("Cannot allocate target while retaining training representation in every stratum.")
  raw <- target * counts_by_stratum / sum(counts_by_stratum)
  allocation <- pmin(floor(raw), capacity)
  while (sum(allocation) < target) {
    available <- which(allocation < capacity)
    if (!length(available)) stop("Allocation failed.")
    residual <- raw[available] - floor(raw[available])
    ord <- available[order(residual, counts_by_stratum[available], decreasing = TRUE)]
    for (idx in ord) {
      if (sum(allocation) >= target) break
      if (allocation[[idx]] < capacity[[idx]]) allocation[[idx]] <- allocation[[idx]] + 1L
    }
  }
  allocation
}

split_stratum <- interaction(tb_metadata$cohort_group, tb_metadata$outcome, drop = TRUE, lex.order = TRUE)
stratum_counts <- table(split_stratum)
TEST_N <- nrow(tb_metadata) - floor(TRAIN_FRACTION * nrow(tb_metadata))
test_allocation <- allocate_integer_target(stratum_counts, TEST_N, retain_at_least_one = TRUE)
final_ice_ids <- character(0)
split_rows <- list()
for (i in seq_along(test_allocation)) {
  stratum_name <- names(test_allocation)[[i]]
  ids <- tb_metadata[split_stratum == stratum_name, sample_id]
  set.seed(seed_for(100, i))
  selected <- if (test_allocation[[i]] > 0L) sample(ids, test_allocation[[i]], replace = FALSE) else character(0)
  final_ice_ids <- c(final_ice_ids, selected)
  split_rows[[i]] <- data.table::data.table(Stratum = stratum_name, Available = length(ids), Final_ice = length(selected))
}
development_ids <- setdiff(tb_metadata$sample_id, final_ice_ids)
development_metadata <- tb_metadata[match(development_ids, sample_id)]
final_ice_metadata <- tb_metadata[match(final_ice_ids, sample_id)]
development_counts <- tb_counts[, development_metadata$sample_id, drop = FALSE]
final_ice_counts <- tb_counts[, final_ice_metadata$sample_id, drop = FALSE]
stopifnot(length(intersect(development_ids, final_ice_ids)) == 0L)
stopifnot(nrow(development_metadata) + nrow(final_ice_metadata) == nrow(tb_metadata))
stopifnot(all(unique(final_ice_metadata$cohort_group) %in% unique(development_metadata$cohort_group)))

split_summary <- data.table::rbindlist(list(
  data.table::data.table(Partition = "Full population", n = nrow(tb_metadata), Controls = sum(tb_metadata$outcome == NEGATIVE_CLASS), Cases = sum(tb_metadata$outcome == POSITIVE_CLASS), Cohorts = uniqueN(tb_metadata$cohort_group)),
  data.table::data.table(Partition = "Development (80%)", n = nrow(development_metadata), Controls = sum(development_metadata$outcome == NEGATIVE_CLASS), Cases = sum(development_metadata$outcome == POSITIVE_CLASS), Cohorts = uniqueN(development_metadata$cohort_group)),
  data.table::data.table(Partition = "Final Ice (20%)", n = nrow(final_ice_metadata), Controls = sum(final_ice_metadata$outcome == NEGATIVE_CLASS), Cases = sum(final_ice_metadata$outcome == POSITIVE_CLASS), Cohorts = uniqueN(final_ice_metadata$cohort_group))
))
data.table::fwrite(split_summary, file.path(split_dir, "CD_full_train_test_split_summary.csv"))
data.table::fwrite(rbindlist(split_rows), file.path(split_dir, "CD_full_train_test_stratum_allocation.csv"))
data.table::fwrite(development_metadata, file.path(split_dir, "CD_development_metadata.csv"))
data.table::fwrite(final_ice_metadata, file.path(split_dir, "CD_FINAL_ICE_metadata.csv"))
saveRDS(development_counts, file.path(split_dir, "CD_development_raw_counts.rds"))
saveRDS(final_ice_counts, file.path(split_dir, "CD_FINAL_ICE_raw_counts.rds"))
print(split_summary)

cohort_names <- sort(unique(development_metadata$cohort_group))

# ======================================================================
# COHORT x OUTCOME-STRATIFIED FOLD MAP
# ======================================================================

make_stratified_fold_map <- function(metadata, folds, repeats, seed) {
  output_rows <- list(); counter <- 0L
  for (repeat_index in seq_len(repeats)) {
    strata <- interaction(metadata$cohort_group, metadata$outcome, drop = TRUE, lex.order = TRUE)
    fold_assignment <- integer(nrow(metadata))
    for (stratum_name in levels(strata)) {
      idx <- which(strata == stratum_name)
      set.seed(seed_for(seed, repeat_index, match(stratum_name, levels(strata))))
      idx <- sample(idx, length(idx), replace = FALSE)
      fold_assignment[idx] <- rep(seq_len(folds), length.out = length(idx))
    }
    for (i in seq_len(nrow(metadata))) {
      counter <- counter + 1L
      output_rows[[counter]] <- data.table::data.table(
        sample_id = metadata$sample_id[[i]], cohort_group = metadata$cohort_group[[i]],
        outcome = as.character(metadata$outcome[[i]]), Repeat = repeat_index, Fold = fold_assignment[[i]]
      )
    }
  }
  rbindlist(output_rows)
}

fold_map <- make_stratified_fold_map(development_metadata, CV_FOLDS, CV_REPEATS, SEED)
data.table::fwrite(fold_map, file.path(benchmark_dir, "CD_development_CV_fold_map.csv"))


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
        "UC_development_",
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
    "CD_development_representation_state_audit.csv"
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
    "CD_representation_diagnostics.csv"
  )
)

data.table::fwrite(
  silhouette_curves,
  file.path(
    diagnostic_dir,
    "CD_representation_silhouette_curves.csv"
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
    "CD_sequencing_depth_dependence.csv"
  )
)

data.table::fwrite(
  technical_by_cohort,
  file.path(
    diagnostic_dir,
    "CD_sequencing_depth_by_cohort.csv"
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
    "CD_feature_sets_full_development.csv"
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
        "CD_Jaccard_matrix_",
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
    "CD_Jaccard_pairwise.csv"
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


# ======================================================================
# FULL-DATA CBA HELPERS
# ======================================================================

compute_class_weights <- function(outcome) {
  y <- factor(as.character(outcome), levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
  counts <- table(y)
  n <- length(y)
  weights <- n / (2 * as.numeric(counts))
  names(weights) <- names(counts)
  weights
}

cohort_stratified_smote <- function(expression_matrix, metadata, seed, k = SMOTE_K) {
  x <- t(expression_matrix)
  y <- as.character(metadata$outcome)
  cohort <- as.character(metadata$cohort_group)
  synthetic_x <- list(); synthetic_meta <- list(); counter <- 0L
  set.seed(seed)
  for (cohort_name in unique(cohort)) {
    idx_c <- which(cohort == cohort_name)
    tab <- table(y[idx_c])
    if (length(tab) < 2L) next
    minority <- names(tab)[which.min(tab)]
    majority_n <- max(tab)
    minority_idx <- idx_c[y[idx_c] == minority]
    need <- majority_n - length(minority_idx)
    if (need <= 0L || length(minority_idx) < 2L) next
    k_use <- min(k, length(minority_idx) - 1L)
    minority_x <- x[minority_idx, , drop = FALSE]
    d <- as.matrix(stats::dist(minority_x))
    diag(d) <- Inf
    neighbours <- lapply(seq_len(nrow(minority_x)), function(i) order(d[i, ])[seq_len(k_use)])
    for (j in seq_len(need)) {
      base_local <- sample(seq_len(nrow(minority_x)), 1L)
      neighbour_local <- sample(neighbours[[base_local]], 1L)
      lambda <- stats::runif(1L)
      new_row <- minority_x[base_local, ] + lambda * (minority_x[neighbour_local, ] - minority_x[base_local, ])
      counter <- counter + 1L
      synthetic_x[[counter]] <- new_row
      synthetic_meta[[counter]] <- data.table::data.table(
        sample_id = paste0("SMOTE_", cohort_name, "_", minority, "_", seed, "_", counter),
        cohort_group = cohort_name,
        outcome = factor(minority, levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)),
        disease = "RA"
      )
    }
  }
  if (counter == 0L) return(list(matrix = expression_matrix, metadata = copy(metadata), synthetic_n = 0L))
  sx <- do.call(rbind, synthetic_x)
  colnames(sx) <- colnames(x)
  augmented_x <- rbind(x, sx)
  augmented_meta <- rbindlist(c(list(copy(metadata)), synthetic_meta), use.names = TRUE, fill = TRUE)
  stopifnot(nrow(augmented_x) == nrow(augmented_meta))
  list(matrix = t(augmented_x), metadata = augmented_meta, synthetic_n = counter)
}

apply_cba_training <- function(expression_matrix, metadata, cba_id, seed) {
  if (cba_id == "None") return(list(matrix = expression_matrix, metadata = copy(metadata), class_weights = NULL, synthetic_n = 0L))
  if (cba_id == "Cohort_SMOTE") {
    sm <- cohort_stratified_smote(expression_matrix, metadata, seed)
    return(list(matrix = sm$matrix, metadata = sm$metadata, class_weights = NULL, synthetic_n = sm$synthetic_n))
  }
  if (cba_id == "Class_Weighting") {
    return(list(matrix = expression_matrix, metadata = copy(metadata), class_weights = compute_class_weights(metadata$outcome), synthetic_n = 0L))
  }
  stop("Unknown CBA: ", cba_id)
}

make_tb_learner <- function(learner_id, class_weights = NULL, hyperparameters = list()) {
  if (learner_id == "RBF_SVM") {
    args <- c(list("classif.svm", predict_type = "prob", type = "C-classification", kernel = "radial", scale = TRUE,
                   cachesize = 500, tolerance = 1e-3, shrinking = TRUE), hyperparameters)
    if (!is.null(class_weights)) args$class.weights <- class_weights
    return(do.call(mlr3::lrn, args))
  }
  if (learner_id == "Linear_SVM") {
    args <- c(list("classif.svm", predict_type = "prob", type = "C-classification", kernel = "linear", scale = TRUE,
                   cachesize = 500, tolerance = 1e-3, shrinking = TRUE), hyperparameters)
    if (!is.null(class_weights)) args$class.weights <- class_weights
    return(do.call(mlr3::lrn, args))
  }
  if (learner_id == "Random_Forest") {
    args <- c(list("classif.ranger", predict_type = "prob", num.trees = 500L, num.threads = 1L,
                   importance = "none", save.memory = TRUE), hyperparameters)
    if (!is.null(class_weights)) args$class.weights <- class_weights
    return(do.call(mlr3::lrn, args))
  }
  stop("Unknown learner: ", learner_id)
}

# ======================================================================
# FAST LEAKAGE-FREE BENCHMARK: 3 x 3 x 3 x 3 = 81 PIPELINES
# ======================================================================

benchmark_rows <- list(); benchmark_counter <- 0L
for (repeat_index in seq_len(CV_REPEATS)) {
  for (fold_index in seq_len(CV_FOLDS)) {
    validation_ids <- fold_map[Repeat == repeat_index & Fold == fold_index, sample_id]
    train_ids <- setdiff(development_metadata$sample_id, validation_ids)
    train_metadata <- development_metadata[match(train_ids, sample_id)]
    validation_metadata <- development_metadata[match(validation_ids, sample_id)]
    train_counts <- development_counts[, train_metadata$sample_id, drop = FALSE]
    validation_counts <- development_counts[, validation_metadata$sample_id, drop = FALSE]
    message("Benchmark fold ", fold_index, "/", CV_FOLDS)
    for (representation_index in seq_along(representation_ids)) {
      representation_id <- representation_ids[[representation_index]]
      rep_fit <- tryCatch(
        fit_apply_representation(representation_id, train_counts, validation_counts, train_metadata, validation_metadata,
                                 seed_for(400, repeat_index, fold_index, representation_index)),
        error = function(e) e
      )
      if (inherits(rep_fit, "error")) next
      for (fsf_index in seq_along(fsf_ids)) {
        fsf_id <- fsf_ids[[fsf_index]]
        fs_fit <- tryCatch(
          select_top_features(rep_fit$train, rep_fit$new, train_metadata, fsf_id, TOP_N,
                              seed_for(500, repeat_index, fold_index, representation_index, fsf_index)),
          error = function(e) e
        )
        if (inherits(fs_fit, "error")) next
        for (cba_index in seq_along(cba_ids)) {
          cba_id <- cba_ids[[cba_index]]
          cba_fit <- apply_cba_training(fs_fit$train, train_metadata, cba_id,
                                        seed_for(550, repeat_index, fold_index, representation_index, fsf_index, cba_index))
          train_task <- make_classification_task(cba_fit$matrix, cba_fit$metadata,
                                                 paste0("BENCH_train_", fold_index, "_", representation_id, "_", fsf_id, "_", cba_id))
          validation_task <- make_classification_task(fs_fit$validation, validation_metadata,
                                                      paste0("BENCH_val_", fold_index, "_", representation_id, "_", fsf_id, "_", cba_id))
          for (learner_index in seq_along(learner_ids)) {
            learner_id <- learner_ids[[learner_index]]
            result <- tryCatch({
              set.seed(seed_for(600, repeat_index, fold_index, representation_index, fsf_index, cba_index, learner_index))
              learner <- make_tb_learner(learner_id, cba_fit$class_weights)
              learner$train(train_task)
              prediction <- learner$predict(validation_task)
              metrics <- safe_prediction_metrics(prediction)
              data.table::data.table(
                Repeat = repeat_index, Fold = fold_index,
                Representation = representation_labels[[representation_id]], Representation_ID = representation_id,
                FS_method = fsf_labels[[fsf_id]], FSF_ID = fsf_id,
                CBA = cba_labels[[cba_id]], CBA_ID = cba_id,
                Learner = learner_labels[[learner_id]], Learner_ID = learner_id,
                Train_n = nrow(train_metadata), Validation_n = nrow(validation_metadata), Synthetic_n = cba_fit$synthetic_n,
                AUROC = metrics$AUROC, AUPRC = metrics$AUPRC, BACC = metrics$BACC, Brier = metrics$Brier,
                Status = "success", Error = NA_character_
              )
            }, error = function(e) data.table::data.table(
              Repeat = repeat_index, Fold = fold_index,
              Representation = representation_labels[[representation_id]], Representation_ID = representation_id,
              FS_method = fsf_labels[[fsf_id]], FSF_ID = fsf_id,
              CBA = cba_labels[[cba_id]], CBA_ID = cba_id,
              Learner = learner_labels[[learner_id]], Learner_ID = learner_id,
              Train_n = nrow(train_metadata), Validation_n = nrow(validation_metadata), Synthetic_n = cba_fit$synthetic_n,
              AUROC = NA_real_, AUPRC = NA_real_, BACC = NA_real_, Brier = NA_real_,
              Status = "model_failed", Error = conditionMessage(e)
            ))
            benchmark_counter <- benchmark_counter + 1L
            benchmark_rows[[benchmark_counter]] <- result
          }
        }
      }
    }
    fwrite(rbindlist(benchmark_rows, use.names = TRUE, fill = TRUE), file.path(benchmark_dir, "CD_full_benchmark_fold_results_PARTIAL.csv"))
  }
}

benchmark_results <- rbindlist(benchmark_rows, use.names = TRUE, fill = TRUE)
fwrite(benchmark_results, file.path(benchmark_dir, "CD_full_benchmark_fold_results.csv"))

benchmark_summary <- benchmark_results[Status == "success" & is.finite(BACC), .(
  Completed_folds = .N,
  Mean_BACC = mean(BACC), SD_BACC = sd(BACC), Worst_BACC = min(BACC),
  Mean_AUROC = mean(AUROC), SD_AUROC = sd(AUROC), Worst_AUROC = min(AUROC),
  Mean_AUPRC = mean(AUPRC), Mean_Brier = mean(Brier), Mean_synthetic_n = mean(Synthetic_n)
), by = .(Representation, Representation_ID, FS_method, FSF_ID, CBA, CBA_ID, Learner, Learner_ID)]
benchmark_summary[, Complete := Completed_folds == CV_FOLDS * CV_REPEATS]
complete_benchmark <- benchmark_summary[Complete == TRUE]
setorder(complete_benchmark, -Mean_BACC, -Worst_BACC, SD_BACC, -Mean_AUROC, Mean_Brier)
complete_benchmark[, Benchmark_rank := seq_len(.N)]
fwrite(complete_benchmark, file.path(benchmark_dir, "CD_full_benchmark_summary.csv"))
if (!nrow(complete_benchmark)) stop("No complete benchmark pipeline.")
print(head(complete_benchmark, 15L))
