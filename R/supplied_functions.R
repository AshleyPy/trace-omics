# Verbatim functions from legacy/equal_size/01.R. Loaded into a private environment.

seed_for <- function(...) {
  values <- as.numeric(c(...))
  weights <- seq_along(values) * 1009
  candidate <- as.numeric(SEED) + sum(values * weights)
  as.integer(candidate %% 2000000000)
}

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
