# ======================================================================
# UC FAST RESCUE MODEL
#
# FAST, LEAKAGE-FREE MODEL SELECTION ON DEVELOPMENT DATA
#
# REPRESENTATIONS:
#   1. Frozen VST only
#   2. log2(CPM + 1)
#
# FEATURE SELECTION:
#   Training-only RF importance
#
# MODELS:
#   RBF SVM
#   Random Forest
#
# PRIMARY SELECTION:
#   Worst strict-LOCO AUROC
#
# FINAL:
#   Fit on all 372 development samples
#   Apply once to the existing 93-sample ICE
#
# IMPORTANT:
#   This creates a NEW directory:
#   TB_FULL_DATA_FINAL_MODEL/101_FAST_RESCUE
# ======================================================================


# ======================================================================
# 0. CLEAN SESSION
# ======================================================================

rm(
  list = ls(
    all.names = TRUE
  )
)

gc(
  verbose = FALSE
)


# ======================================================================
# 1. PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "DESeq2",
  "SummarizedExperiment",
  "matrixStats",
  "ranger",
  "e1071",
  "pROC"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (
  length(
    missing_packages
  ) > 0L
) {

  stop(
    "Install these packages first:\n",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({

  library(
    data.table
  )
})


# ======================================================================
# 2. SETTINGS
# ======================================================================

SEED <- 7832L

POSITIVE_CLASS <- "case"
NEGATIVE_CLASS <- "control"

TOP_N <- 50L

MIN_TOTAL_TRAINING_COUNT <- 10L

# Faster than running importance across every gene.
# This is training-only and therefore leakage-free.
RF_PREFILTER <- 3000L

RF_IMPORTANCE_TREES <- 250L

RF_MODEL_TREES <- 1000L

MIN_ACCEPTABLE_WORST_AUROC <- 0.65


PROJECT_DIR <- paste0(
  Sys.getenv(
    "USERPROFILE"
  ),
  "/Desktop/UC_FULL_DATA_FINAL_MODEL"
)

SPLIT_DIR <- file.path(
  PROJECT_DIR,
  "01_split"
)

RESCUE_DIR <- file.path(
  PROJECT_DIR,
  "101_FAST_RESCUE"
)

dir.create(
  RESCUE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 3. LOAD EXISTING DEVELOPMENT / ICE SPLIT
# ======================================================================

development_counts <- readRDS(
  file.path(
    SPLIT_DIR,
    "UC_development_raw_counts.rds"
  )
)

final_ice_counts <- readRDS(
  file.path(
    SPLIT_DIR,
    "UC_FINAL_ICE_raw_counts.rds"
  )
)

development_metadata <- data.table::fread(
  file.path(
    SPLIT_DIR,
    "UC_development_metadata.csv"
  )
)

final_ice_metadata <- data.table::fread(
  file.path(
    SPLIT_DIR,
    "UC_FINAL_ICE_metadata.csv"
  )
)


# ======================================================================
# 4. STANDARDISE METADATA
# ======================================================================

prepare_metadata <- function(
  metadata
) {

  metadata <- data.table::as.data.table(
    metadata
  )

  if (
    !"sample_id" %in%
      names(
        metadata
      )
  ) {

    metadata[
      ,
      sample_id :=
        as.character(
          geo_accession
        )
    ]
  }

  if (
    !"cohort_group" %in%
      names(
        metadata
      )
  ) {

    metadata[
      ,
      cohort_group :=
        as.character(
          cohort
        )
    ]
  }

  metadata[
    ,
    sample_id :=
      trimws(
        as.character(
          sample_id
        )
      )
  ]

  metadata[
    ,
    cohort_group :=
      trimws(
        as.character(
          cohort_group
        )
      )
  ]

  metadata[
    ,
    outcome :=
      factor(
        tolower(
          trimws(
            as.character(
              outcome
            )
          )
        ),
        levels = c(
          NEGATIVE_CLASS,
          POSITIVE_CLASS
        )
      )
  ]

  if (
    anyNA(
      metadata$outcome
    )
  ) {

    stop(
      "Invalid outcome label detected."
    )
  }

  metadata
}


development_metadata <- prepare_metadata(
  development_metadata
)

final_ice_metadata <- prepare_metadata(
  final_ice_metadata
)


# ======================================================================
# 5. ALIGN MATRICES
# ======================================================================

development_counts <- as.matrix(
  development_counts
)

final_ice_counts <- as.matrix(
  final_ice_counts
)

storage.mode(
  development_counts
) <- "integer"

storage.mode(
  final_ice_counts
) <- "integer"


if (
  !identical(
    rownames(
      development_counts
    ),
    rownames(
      final_ice_counts
    )
  )
) {

  stop(
    "Development and ICE gene order differs."
  )
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


stopifnot(
  identical(
    colnames(
      development_counts
    ),
    development_metadata$sample_id
  )
)

stopifnot(
  identical(
    colnames(
      final_ice_counts
    ),
    final_ice_metadata$sample_id
  )
)


cat(
  "\n============================================================\n",
  "DATA READY\n",
  "============================================================\n",
  "Development: ",
  nrow(
    development_metadata
  ),
  "\n",
  "ICE: ",
  nrow(
    final_ice_metadata
  ),
  "\n",
  "Genes: ",
  nrow(
    development_counts
  ),
  "\n",
  "Cases: ",
  sum(
    development_metadata$outcome ==
      POSITIVE_CLASS
  ),
  "\n",
  "Controls: ",
  sum(
    development_metadata$outcome ==
      NEGATIVE_CLASS
  ),
  "\n",
  sep = ""
)


# ======================================================================
# 6. SEED HELPER
# ======================================================================

seed_for <- function(...) {

  values <- as.numeric(
    c(...)
  )

  weights <- seq_along(
    values
  ) * 1009

  as.integer(
    (
      SEED +
        sum(
          values *
            weights
        )
    ) %%
      2000000000
  )
}


# ======================================================================
# 7. CLASS WEIGHTS
# ======================================================================

compute_class_weights <- function(
  outcome
) {

  y <- factor(
    as.character(
      outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  counts <- table(
    y
  )

  n <- length(
    y
  )

  weights <- n /
    (
      2 *
        as.numeric(
          counts
        )
    )

  names(
    weights
  ) <- names(
    counts
  )

  weights
}


# ======================================================================
# 8. METRICS
# ======================================================================

calculate_metrics <- function(
  truth,
  probability,
  threshold = 0.5
) {

  truth <- factor(
    as.character(
      truth
    ),
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

  response <- ifelse(
    probability >=
      threshold,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  )

  sensitivity <- mean(
    response[
      truth ==
        POSITIVE_CLASS
    ] ==
      POSITIVE_CLASS
  )

  specificity <- mean(
    response[
      truth ==
        NEGATIVE_CLASS
    ] ==
      NEGATIVE_CLASS
  )

  bacc <- mean(
    c(
      sensitivity,
      specificity
    )
  )

  numeric_truth <- as.integer(
    truth ==
      POSITIVE_CLASS
  )

  list(
    AUROC =
      as.numeric(
        pROC::auc(
          roc_object
        )
      ),
    BACC =
      bacc,
    Sensitivity =
      sensitivity,
    Specificity =
      specificity,
    Brier =
      mean(
        (
          probability -
            numeric_truth
        )^2
      )
  )
}


calculate_bacc <- function(
  truth,
  probability,
  threshold
) {

  truth <- factor(
    as.character(
      truth
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  response <- ifelse(
    probability >=
      threshold,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  )

  sensitivity <- mean(
    response[
      truth ==
        POSITIVE_CLASS
    ] ==
      POSITIVE_CLASS
  )

  specificity <- mean(
    response[
      truth ==
        NEGATIVE_CLASS
    ] ==
      NEGATIVE_CLASS
  )

  mean(
    c(
      sensitivity,
      specificity
    )
  )
}


# ======================================================================
# 9. FROZEN VST ONLY
#
# This retains your training-derived VST strategy but removes fSVA.
# ======================================================================

calculate_poscounts_geo_means <- function(
  count_matrix
) {

  positive_number <- rowSums(
    count_matrix >
      0
  )

  if (
    any(
      positive_number ==
        0L
    )
  ) {

    stop(
      "Zero-only training gene detected."
    )
  }

  log_counts <- log(
    count_matrix
  )

  log_counts[
    !is.finite(
      log_counts
    )
  ] <- 0

  geo_means <- exp(
    rowSums(
      log_counts
    ) /
      positive_number
  )

  names(
    geo_means
  ) <- rownames(
    count_matrix
  )

  geo_means
}


fit_apply_vst_only <- function(
  train_counts,
  new_counts = NULL,
  train_metadata
) {

  train_counts <- train_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]


  if (
    !is.null(
      new_counts
    )
  ) {

    new_counts <- new_counts[
      rownames(
        train_counts
      ),
      ,
      drop = FALSE
    ]
  }


  training_total_counts <- rowSums(
    train_counts
  )


  keep_gene <- (
    is.finite(
      training_total_counts
    ) &
      training_total_counts >=
        MIN_TOTAL_TRAINING_COUNT
  )


  train_filtered <- train_counts[
    keep_gene,
    ,
    drop = FALSE
  ]


  new_filtered <- NULL


  if (
    !is.null(
      new_counts
    )
  ) {

    new_filtered <- new_counts[
      keep_gene,
      ,
      drop = FALSE
    ]
  }


  train_coldata <- data.frame(
    outcome =
      factor(
        as.character(
          train_metadata$outcome
        ),
        levels = c(
          NEGATIVE_CLASS,
          POSITIVE_CLASS
        )
      ),
    row.names =
      train_metadata$sample_id
  )


  geo_means <- calculate_poscounts_geo_means(
    train_filtered
  )


  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData =
      train_filtered,
    colData =
      train_coldata,
    design =
      ~ outcome
  )


  dds_train <- DESeq2::estimateSizeFactors(
    dds_train,
    type = "poscounts",
    geoMeans =
      geo_means
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


  dispersion_function <- DESeq2::dispersionFunction(
    dds_train
  )


  vst_train <- DESeq2::varianceStabilizingTransformation(
    dds_train,
    blind = FALSE
  )


  train_matrix <- as.matrix(
    SummarizedExperiment::assay(
      vst_train
    )
  )


  new_matrix <- NULL


  if (
    !is.null(
      new_filtered
    )
  ) {

    new_ids <- colnames(
      new_filtered
    )


    new_coldata <- data.frame(
      intercept =
        rep(
          1L,
          length(
            new_ids
          )
        ),
      row.names =
        new_ids
    )


    dds_new <- DESeq2::DESeqDataSetFromMatrix(
      countData =
        new_filtered,
      colData =
        new_coldata,
      design =
        ~ 1
    )


    dds_new <- DESeq2::estimateSizeFactors(
      dds_new,
      type = "poscounts",
      geoMeans =
        geo_means
    )


    DESeq2::dispersionFunction(
      dds_new
    ) <- dispersion_function


    vst_new <- DESeq2::varianceStabilizingTransformation(
      dds_new,
      blind = FALSE
    )


    new_matrix <- as.matrix(
      SummarizedExperiment::assay(
        vst_new
      )
    )
  }


  list(
    train =
      train_matrix,
    new =
      new_matrix,
    state =
      list(
        raw_genes =
          rownames(
            train_counts
          ),
        genes =
          rownames(
            train_filtered
          ),
        geo_means =
          geo_means,
        dispersion_function =
          dispersion_function
      )
  )
}


# ======================================================================
# 10. SAMPLE-WISE log2(CPM + 1)
#
# Very simple deployment path.
# No fSVA projection.
# ======================================================================

logcpm_transform <- function(
  counts
) {

  library_size <- colSums(
    counts
  )


  if (
    any(
      library_size <=
        0
    )
  ) {

    stop(
      "Invalid RNA-seq library size."
    )
  }


  cpm <- sweep(
    counts,
    MARGIN = 2L,
    STATS =
      library_size,
    FUN = "/"
  ) *
    1e6


  log2(
    cpm +
      1
  )
}


fit_apply_logcpm <- function(
  train_counts,
  new_counts = NULL,
  train_metadata
) {

  train_counts <- train_counts[
    ,
    train_metadata$sample_id,
    drop = FALSE
  ]


  if (
    !is.null(
      new_counts
    )
  ) {

    new_counts <- new_counts[
      rownames(
        train_counts
      ),
      ,
      drop = FALSE
    ]
  }


  training_total_counts <- rowSums(
    train_counts
  )


  keep_gene <- (
    is.finite(
      training_total_counts
    ) &
      training_total_counts >=
        MIN_TOTAL_TRAINING_COUNT
  )


  train_matrix <- logcpm_transform(
    train_counts
  )[
    keep_gene,
    ,
    drop = FALSE
  ]


  new_matrix <- NULL


  if (
    !is.null(
      new_counts
    )
  ) {

    new_matrix <- logcpm_transform(
      new_counts
    )[
      keep_gene,
      ,
      drop = FALSE
    ]
  }


  list(
    train =
      train_matrix,
    new =
      new_matrix,
    state =
      list(
        raw_genes =
          rownames(
            train_counts
          ),
        genes =
          rownames(
            train_counts
          )[
            keep_gene
          ]
      )
  )
}


# ======================================================================
# 11. REPRESENTATION WRAPPER
# ======================================================================

fit_apply_representation <- function(
  representation,
  train_counts,
  new_counts,
  train_metadata
) {

  if (
    representation ==
      "VST_only"
  ) {

    return(
      fit_apply_vst_only(
        train_counts =
          train_counts,
        new_counts =
          new_counts,
        train_metadata =
          train_metadata
      )
    )
  }


  if (
    representation ==
      "logCPM"
  ) {

    return(
      fit_apply_logcpm(
        train_counts =
          train_counts,
        new_counts =
          new_counts,
        train_metadata =
          train_metadata
      )
    )
  }


  stop(
    "Unknown representation."
  )
}


# ======================================================================
# 12. FAST TRAINING-ONLY RF IMPORTANCE
#
# First restrict to top variance genes WITHIN TRAINING DATA ONLY.
# Then run RF permutation importance.
# ======================================================================

select_rf_features <- function(
  train_matrix,
  train_metadata,
  seed
) {

  gene_variance <- matrixStats::rowVars(
    train_matrix
  )


  valid <- (
    is.finite(
      gene_variance
    ) &
      gene_variance >
        0
  )


  train_matrix <- train_matrix[
    valid,
    ,
    drop = FALSE
  ]


  gene_variance <- gene_variance[
    valid
  ]


  if (
    nrow(
      train_matrix
    ) >
      RF_PREFILTER
  ) {

    variance_order <- order(
      gene_variance,
      decreasing = TRUE
    )


    keep <- variance_order[
      seq_len(
        RF_PREFILTER
      )
    ]


    train_matrix <- train_matrix[
      keep,
      ,
      drop = FALSE
    ]
  }


  x <- as.data.frame(
    t(
      train_matrix
    ),
    check.names = FALSE
  )


  y <- factor(
    as.character(
      train_metadata$outcome
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )


  class_weights <- compute_class_weights(
    y
  )


  set.seed(
    seed
  )


  rf <- ranger::ranger(
    x = x,
    y = y,
    probability = TRUE,
    importance = "permutation",
    num.trees =
      RF_IMPORTANCE_TREES,
    mtry =
      max(
        1L,
        floor(
          sqrt(
            ncol(
              x
            )
          )
        )
      ),
    class.weights =
      class_weights,
    num.threads = 0L,
    seed = seed
  )


  importance_values <- rf$variable.importance


  importance_values[
    !is.finite(
      importance_values
    )
  ] <- -Inf


  selected_genes <- names(
    head(
      sort(
        importance_values,
        decreasing = TRUE
      ),
      TOP_N
    )
  )


  if (
    length(
      selected_genes
    ) !=
      TOP_N
  ) {

    stop(
      "RF feature selection did not return 50 genes."
    )
  }


  selected_genes
}


# ======================================================================
# 13. MODEL CANDIDATES
#
# Only four useful RBF settings + one Random Forest.
# This is deliberately small for speed.
# ======================================================================

candidate_table <- data.table::rbindlist(
  list(

    data.table::data.table(
      Model = "RBF SVM",
      Config = "RBF_1",
      cost = 1,
      gamma = 0.0200000
    ),

    data.table::data.table(
      Model = "RBF SVM",
      Config = "RBF_2",
      cost = 3.5345224,
      gamma = 0.001012607
    ),

    data.table::data.table(
      Model = "RBF SVM",
      Config = "RBF_3",
      cost = 1.6746829,
      gamma = 0.0003574437
    ),

    data.table::data.table(
      Model = "RBF SVM",
      Config = "RBF_4",
      cost = 0.5606065,
      gamma = 0.007523222
    ),

    data.table::data.table(
      Model = "Random Forest",
      Config = "RF_1",
      cost = NA_real_,
      gamma = NA_real_
    )
  )
)


# ======================================================================
# 14. FIT MODEL
# ======================================================================

fit_model <- function(
  x_train,
  y_train,
  candidate,
  seed
) {

  y_train <- factor(
    as.character(
      y_train
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )


  class_weights <- compute_class_weights(
    y_train
  )


  if (
    candidate$Model[[1L]] ==
      "RBF SVM"
  ) {

    set.seed(
      seed
    )


    model <- e1071::svm(
      x = x_train,
      y = y_train,
      type = "C-classification",
      kernel = "radial",
      cost =
        candidate$cost[[1L]],
      gamma =
        candidate$gamma[[1L]],
      class.weights =
        class_weights,
      probability = TRUE,
      scale = TRUE,
      cachesize = 500,
      tolerance = 0.001,
      shrinking = TRUE
    )


    return(
      list(
        type = "RBF SVM",
        fit = model
      )
    )
  }


  if (
    candidate$Model[[1L]] ==
      "Random Forest"
  ) {

    set.seed(
      seed
    )


    model <- ranger::ranger(
      x =
        as.data.frame(
          x_train,
          check.names = FALSE
        ),
      y =
        y_train,
      probability = TRUE,
      num.trees =
        RF_MODEL_TREES,
      mtry =
        max(
          1L,
          floor(
            sqrt(
              ncol(
                x_train
              )
            )
          )
        ),
      min.node.size = 5L,
      class.weights =
        class_weights,
      importance = "none",
      num.threads = 0L,
      seed = seed
    )


    return(
      list(
        type = "Random Forest",
        fit = model
      )
    )
  }


  stop(
    "Unknown model."
  )
}


# ======================================================================
# 15. PREDICT MODEL
# ======================================================================

predict_model <- function(
  model,
  x_test
) {

  if (
    model$type ==
      "RBF SVM"
  ) {

    prediction <- stats::predict(
      model$fit,
      newdata =
        x_test,
      probability = TRUE
    )


    probability_matrix <- attr(
      prediction,
      "probabilities"
    )


    return(
      as.numeric(
        probability_matrix[
          ,
          POSITIVE_CLASS
        ]
      )
    )
  }


  if (
    model$type ==
      "Random Forest"
  ) {

    probability_matrix <- predict(
      model$fit,
      data =
        as.data.frame(
          x_test,
          check.names = FALSE
        )
    )$predictions


    return(
      as.numeric(
        probability_matrix[
          ,
          POSITIVE_CLASS
        ]
      )
    )
  }


  stop(
    "Unknown fitted model."
  )
}


# ======================================================================
# 16. STRICT LOCO MODEL SELECTION
#
# NO ICE DATA ARE USED HERE.
# ======================================================================

representations <- c(
  "VST_only",
  "logCPM"
)


cohort_names <- sort(
  unique(
    development_metadata$cohort_group
  )
)


metric_rows <- list()

prediction_rows <- list()

metric_counter <- 0L

prediction_counter <- 0L


for (
  representation_index in seq_along(
    representations
  )
) {

  representation <- representations[[representation_index]]


  for (
    target_index in seq_along(
      cohort_names
    )
  ) {

    target_cohort <- cohort_names[[target_index]]


    cat(
      "\n============================================================\n",
      representation,
      " | HELD OUT: ",
      target_cohort,
      "\n",
      "============================================================\n",
      sep = ""
    )


    source_metadata <- development_metadata[
      cohort_group !=
        target_cohort
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


    representation_fit <- tryCatch(

      fit_apply_representation(
        representation =
          representation,
        train_counts =
          source_counts,
        new_counts =
          target_counts,
        train_metadata =
          source_metadata
      ),

      error = function(e) {

        message(
          "Representation failed: ",
          conditionMessage(
            e
          )
        )

        NULL
      }
    )


    if (
      is.null(
        representation_fit
      )
    ) {

      next
    }


    selected_genes <- select_rf_features(
      train_matrix =
        representation_fit$train,
      train_metadata =
        source_metadata,
      seed =
        seed_for(
          100,
          representation_index,
          target_index
        )
    )


    x_train <- t(
      representation_fit$train[
        selected_genes,
        ,
        drop = FALSE
      ]
    )


    x_test <- t(
      representation_fit$new[
        selected_genes,
        ,
        drop = FALSE
      ]
    )


    y_train <- source_metadata$outcome

    y_test <- target_metadata$outcome


    for (
      candidate_index in seq_len(
        nrow(
          candidate_table
        )
      )
    ) {

      candidate <- candidate_table[
        candidate_index
      ]


      message(
        "  ",
        candidate$Config[[1L]]
      )


      result <- tryCatch(
        {

          current_model <- fit_model(
            x_train =
              x_train,
            y_train =
              y_train,
            candidate =
              candidate,
            seed =
              seed_for(
                200,
                representation_index,
                target_index,
                candidate_index
              )
          )


          probability <- predict_model(
            model =
              current_model,
            x_test =
              x_test
          )


          metrics <- calculate_metrics(
            truth =
              y_test,
            probability =
              probability,
            threshold =
              0.5
          )


          list(
            probability =
              probability,
            metrics =
              metrics
          )
        },

        error = function(e) {

          message(
            "    FAILED: ",
            conditionMessage(
              e
            )
          )

          NULL
        }
      )


      if (
        is.null(
          result
        )
      ) {

        next
      }


      metric_counter <- metric_counter +
        1L


      metric_rows[[metric_counter]] <- data.table::data.table(
        Representation =
          representation,
        Target_cohort =
          target_cohort,
        Model =
          candidate$Model[[1L]],
        Config =
          candidate$Config[[1L]],
        cost =
          candidate$cost[[1L]],
        gamma =
          candidate$gamma[[1L]],
        AUROC =
          result$metrics$AUROC,
        BACC =
          result$metrics$BACC,
        Sensitivity =
          result$metrics$Sensitivity,
        Specificity =
          result$metrics$Specificity,
        Brier =
          result$metrics$Brier
      )


      prediction_counter <- prediction_counter +
        1L


      prediction_rows[[prediction_counter]] <- data.table::data.table(
        Representation =
          representation,
        Target_cohort =
          target_cohort,
        Model =
          candidate$Model[[1L]],
        Config =
          candidate$Config[[1L]],
        sample_id =
          target_metadata$sample_id,
        truth =
          as.character(
            y_test
          ),
        probability_case =
          result$probability
      )
    }
  }
}


# ======================================================================
# 17. COMBINE LOCO RESULTS
# ======================================================================

loco_metrics <- data.table::rbindlist(
  metric_rows,
  use.names = TRUE,
  fill = TRUE
)


loco_predictions <- data.table::rbindlist(
  prediction_rows,
  use.names = TRUE,
  fill = TRUE
)


data.table::fwrite(
  loco_metrics,
  file.path(
    RESCUE_DIR,
    "01_LOCO_metrics.csv"
  )
)


data.table::fwrite(
  loco_predictions,
  file.path(
    RESCUE_DIR,
    "02_LOCO_predictions.csv"
  )
)


# ======================================================================
# 18. RANK BY TRANSPORTABILITY
# ======================================================================

loco_summary <- loco_metrics[
  ,
  .(
    Cohorts =
      data.table::uniqueN(
        Target_cohort
      ),
    Mean_LOCO_AUROC =
      mean(
        AUROC
      ),
    Worst_cohort_AUROC =
      min(
        AUROC
      ),
    SD_LOCO_AUROC =
      stats::sd(
        AUROC
      ),
    Mean_LOCO_BACC =
      mean(
        BACC
      ),
    Worst_cohort_BACC =
      min(
        BACC
      ),
    Mean_Brier =
      mean(
        Brier
      )
  ),
  by = .(
    Representation,
    Model,
    Config,
    cost,
    gamma
  )
]


loco_summary <- loco_summary[
  Cohorts ==
    length(
      cohort_names
    )
]


data.table::setorder(
  loco_summary,
  -Worst_cohort_AUROC,
  -Mean_LOCO_AUROC,
  SD_LOCO_AUROC,
  -Mean_LOCO_BACC
)


loco_summary[
  ,
  Rank :=
    seq_len(
      .N
    )
]


data.table::fwrite(
  loco_summary,
  file.path(
    RESCUE_DIR,
    "03_LOCO_ranked_summary.csv"
  )
)


cat(
  "\n============================================================\n",
  "STRICT LOCO RANKING\n",
  "============================================================\n"
)


print(
  loco_summary
)


if (
  nrow(
    loco_summary
  ) ==
    0L
) {

  stop(
    "No model completed all three LOCO cohorts."
  )
}


best_model <- loco_summary[
  1L
]


cat(
  "\n============================================================\n",
  "WINNING RESCUE MODEL\n",
  "============================================================\n"
)


print(
  best_model
)


if (
  best_model$Worst_cohort_AUROC[[1L]] <
    MIN_ACCEPTABLE_WORST_AUROC
) {

  warning(
    "Worst-cohort AUROC is below ",
    MIN_ACCEPTABLE_WORST_AUROC,
    ". Do not call this deployable yet."
  )
}


# ======================================================================
# 19. RECOVER SELECTED CANDIDATE
# ======================================================================

selected_representation <- best_model$Representation[[1L]]

selected_model_name <- best_model$Model[[1L]]

selected_config <- best_model$Config[[1L]]


selected_candidate <- candidate_table[
  Config ==
    selected_config
]


# ======================================================================
# 20. TUNE THRESHOLD FROM DEVELOPMENT LOCO PREDICTIONS
#
# This optimises mean BACC across the three unseen cohorts.
# ======================================================================

best_predictions <- loco_predictions[
  Representation ==
    selected_representation &
    Model ==
      selected_model_name &
    Config ==
      selected_config
]


threshold_grid <- seq(
  0.05,
  0.95,
  by = 0.01
)


threshold_rows <- list()


for (
  threshold_index in seq_along(
    threshold_grid
  )
) {

  current_threshold <- threshold_grid[[threshold_index]]


  cohort_bacc <- best_predictions[
    ,
    .(
      BACC =
        calculate_bacc(
          truth =
            truth,
          probability =
            probability_case,
          threshold =
            current_threshold
        )
    ),
    by =
      Target_cohort
  ]


  threshold_rows[[threshold_index]] <- data.table::data.table(
    Threshold =
      current_threshold,
    Mean_cohort_BACC =
      mean(
        cohort_bacc$BACC
      ),
    Worst_cohort_BACC =
      min(
        cohort_bacc$BACC
      )
  )
}


threshold_summary <- data.table::rbindlist(
  threshold_rows
)


threshold_summary[
  ,
  Distance_from_half :=
    abs(
      Threshold -
        0.5
    )
]


data.table::setorder(
  threshold_summary,
  -Mean_cohort_BACC,
  -Worst_cohort_BACC,
  Distance_from_half
)


locked_threshold <- threshold_summary$Threshold[[1L]]


data.table::fwrite(
  threshold_summary,
  file.path(
    RESCUE_DIR,
    "04_threshold_selection.csv"
  )
)


cat(
  "\nLOCKED THRESHOLD: ",
  sprintf(
    "%.2f",
    locked_threshold
  ),
  "\n",
  sep = ""
)


# ======================================================================
# 21. FIT WINNING MODEL ON ALL 372 DEVELOPMENT SAMPLES
# ======================================================================

cat(
  "\n============================================================\n",
  "FINAL DEVELOPMENT FIT\n",
  "============================================================\n"
)


final_representation <- fit_apply_representation(
  representation =
    selected_representation,
  train_counts =
    development_counts,
  new_counts =
    final_ice_counts,
  train_metadata =
    development_metadata
)


final_genes <- select_rf_features(
  train_matrix =
    final_representation$train,
  train_metadata =
    development_metadata,
  seed =
    seed_for(
      900,
      1
    )
)


final_x_train <- t(
  final_representation$train[
    final_genes,
    ,
    drop = FALSE
  ]
)


final_x_ice <- t(
  final_representation$new[
    final_genes,
    ,
    drop = FALSE
  ]
)


final_model <- fit_model(
  x_train =
    final_x_train,
  y_train =
    development_metadata$outcome,
  candidate =
    selected_candidate,
  seed =
    seed_for(
      910,
      1
    )
)


# ======================================================================
# 22. DIAGNOSTIC ICE PREDICTION
# ======================================================================

ice_probability <- predict_model(
  model =
    final_model,
  x_test =
    final_x_ice
)


ice_metrics <- calculate_metrics(
  truth =
    final_ice_metadata$outcome,
  probability =
    ice_probability,
  threshold =
    locked_threshold
)


ice_roc <- pROC::roc(
  response =
    final_ice_metadata$outcome,
  predictor =
    ice_probability,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


ice_ci <- pROC::ci.auc(
  ice_roc,
  method = "delong",
  conf.level = 0.95
)


ice_response <- ifelse(
  ice_probability >=
    locked_threshold,
  POSITIVE_CLASS,
  NEGATIVE_CLASS
)


# ======================================================================
# 23. SAVE RESULTS
# ======================================================================

ice_result <- data.table::data.table(
  Representation =
    selected_representation,
  Model =
    selected_model_name,
  Config =
    selected_config,
  cost =
    selected_candidate$cost[[1L]],
  gamma =
    selected_candidate$gamma[[1L]],
  Threshold =
    locked_threshold,
  Development_n =
    nrow(
      development_metadata
    ),
  ICE_n =
    nrow(
      final_ice_metadata
    ),
  AUROC =
    ice_metrics$AUROC,
  AUROC_lower95 =
    as.numeric(
      ice_ci[[1L]]
    ),
  AUROC_upper95 =
    as.numeric(
      ice_ci[[3L]]
    ),
  BACC =
    ice_metrics$BACC,
  Sensitivity =
    ice_metrics$Sensitivity,
  Specificity =
    ice_metrics$Specificity,
  Brier =
    ice_metrics$Brier
)


ice_predictions <- data.table::data.table(
  sample_id =
    final_ice_metadata$sample_id,
  cohort_group =
    final_ice_metadata$cohort_group,
  truth =
    as.character(
      final_ice_metadata$outcome
    ),
  probability_case =
    ice_probability,
  response =
    ice_response
)


gene_table <- data.table::data.table(
  Rank =
    seq_along(
      final_genes
    ),
  Gene =
    final_genes
)


data.table::fwrite(
  ice_result,
  file.path(
    RESCUE_DIR,
    "05_RESCUE_ICE_result.csv"
  )
)


data.table::fwrite(
  ice_predictions,
  file.path(
    RESCUE_DIR,
    "06_RESCUE_ICE_predictions.csv"
  )
)


data.table::fwrite(
  gene_table,
  file.path(
    RESCUE_DIR,
    "07_RESCUE_final_genes.csv"
  )
)


# ======================================================================
# 24. SAVE WORKING MODEL BUNDLE
# ======================================================================

rescue_bundle <- list(

  bundle_version =
    "UC_FAST_RESCUE_1.0",

  created =
    as.character(
      Sys.time()
    ),

  disease =
    "UC",

  representation =
    selected_representation,

  model_name =
    selected_model_name,

  config =
    selected_config,

  cost =
    selected_candidate$cost[[1L]],

  gamma =
    selected_candidate$gamma[[1L]],

  selected_genes =
    final_genes,

  threshold =
    locked_threshold,

  raw_gene_order =
    rownames(
      development_counts
    ),

  representation_state =
    final_representation$state,

  fitted_model =
    final_model,

  development_sample_ids =
    development_metadata$sample_id,

  strict_loco_selection =
    best_model,

  diagnostic_ice_result =
    ice_result
)


saveRDS(
  rescue_bundle,
  file.path(
    RESCUE_DIR,
    "UC_FAST_RESCUE_MODEL_BUNDLE.rds"
  ),
  compress = "xz"
)


# ======================================================================
# 25. FINAL OUTPUT
# ======================================================================

cat(
  "\n\n============================================================\n",
  "UC FAST RESCUE COMPLETE\n",
  "============================================================\n",
  "Representation: ",
  selected_representation,
  "\n",
  "Model: ",
  selected_model_name,
  "\n",
  "Config: ",
  selected_config,
  "\n",
  "Worst LOCO AUROC: ",
  sprintf(
    "%.3f",
    best_model$Worst_cohort_AUROC[[1L]]
  ),
  "\n",
  "Mean LOCO AUROC: ",
  sprintf(
    "%.3f",
    best_model$Mean_LOCO_AUROC[[1L]]
  ),
  "\n",
  "Threshold: ",
  sprintf(
    "%.2f",
    locked_threshold
  ),
  "\n",
  "\nDIAGNOSTIC ICE\n",
  "AUROC: ",
  sprintf(
    "%.3f",
    ice_result$AUROC[[1L]]
  ),
  " (",
  sprintf(
    "%.3f",
    ice_result$AUROC_lower95[[1L]]
  ),
  "--",
  sprintf(
    "%.3f",
    ice_result$AUROC_upper95[[1L]]
  ),
  ")\n",
  "BACC: ",
  sprintf(
    "%.3f",
    ice_result$BACC[[1L]]
  ),
  "\n",
  "Sensitivity: ",
  sprintf(
    "%.3f",
    ice_result$Sensitivity[[1L]]
  ),
  "\n",
  "Specificity: ",
  sprintf(
    "%.3f",
    ice_result$Specificity[[1L]]
  ),
  "\n",
  "\nSaved to:\n",
  RESCUE_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)
