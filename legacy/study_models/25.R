# ======================================================================
# PART 2
# LIGHTGBM vs TabPFN-3
#
# POOLED NESTED-CV BENCHMARK
#
# DATA:
#   All 372 development samples
#
# NO LOCO.
#
# OUTER:
#   5-fold CV
#   stratified by cohort x outcome
#
# INNER:
#   3-fold CV
#   LightGBM tuning only
#
# REPRESENTATION:
#   log2(CPM + 1)
#
# FEATURE SELECTION:
#   RF permutation importance
#   TOP 100
#   performed independently inside resampling
#
# MODELS:
#   1. tuned LightGBM
#   2. default TabPFN-3 foundation model
#
# PRIMARY METRIC:
#   AUROC
# ======================================================================


# ======================================================================
# 2.1 PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "pROC",
  "PRROC",
  "reticulate"
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
  install.packages(
    missing_packages
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
  library(reticulate)
})


# ======================================================================
# 2.2 AUDIT REQUIRED OBJECTS
# ======================================================================

required_objects <- c(
  "development_counts",
  "development_metadata",
  "fit_apply_logcpm",
  "select_rf_features_final",
  "POSITIVE_CLASS",
  "NEGATIVE_CLASS",
  "PROJECT_DIR"
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
    "Missing required objects:\n",
    paste(
      missing_objects,
      collapse = "\n"
    )
  )
}


# ======================================================================
# 2.3 SETTINGS
# ======================================================================

SEED_FOUNDATION <- 7832L

OUTER_FOLDS <- 5L

INNER_FOLDS <- 3L

FOUNDATION_TOP_N <- 100L


FOUNDATION_DIR <- file.path(
  PROJECT_DIR,
  "103_LIGHTGBM_VS_TABPFN"
)

dir.create(
  FOUNDATION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ======================================================================
# PART 3
# PYTHON MODEL BACKEND
# ======================================================================

TB_PYTHON_ENV <- "tb_foundation_models"


reticulate::use_virtualenv(
  TB_PYTHON_ENV,
  required = TRUE
)


if (
  !reticulate::py_module_available(
    "lightgbm"
  )
) {

  stop(
    "Python package 'lightgbm' is not installed."
  )
}


if (
  !reticulate::py_module_available(
    "tabpfn"
  )
) {

  stop(
    "Python package 'tabpfn' is not installed."
  )
}


np <- reticulate::import(
  "numpy",
  convert = FALSE
)

lightgbm_py <- reticulate::import(
  "lightgbm",
  convert = FALSE
)

tabpfn_py <- reticulate::import(
  "tabpfn",
  convert = FALSE
)

torch_py <- reticulate::import(
  "torch",
  convert = TRUE
)


TABPFN_DEVICE <- if (
  isTRUE(
    torch_py$cuda$is_available()
  )
) {
  "cuda"
} else {
  "cpu"
}


cat(
  "\nTabPFN device: ",
  TABPFN_DEVICE,
  "\n",
  sep = ""
)

# ======================================================================
# PART 4
# STRATIFIED FOLD GENERATOR
#
# Stratifies jointly by:
#
#   cohort x outcome
#
# This keeps all three cohorts represented in each fold.
# ======================================================================

make_joint_stratified_folds <- function(
  metadata,
  number_of_folds,
  seed
) {

  metadata <- data.table::copy(
    metadata
  )

  metadata[
    ,
    stratum := paste(
      cohort_group,
      outcome,
      sep = "__"
    )
  ]

  fold_id <- integer(
    nrow(
      metadata
    )
  )

  set.seed(
    seed
  )

  for (
    current_stratum in unique(
      metadata$stratum
    )
  ) {

    current_index <- which(
      metadata$stratum ==
        current_stratum
    )

    current_index <- sample(
      current_index
    )

    assigned_folds <- rep(
      seq_len(
        number_of_folds
      ),
      length.out = length(
        current_index
      )
    )

    fold_id[
      current_index
    ] <- assigned_folds
  }

  fold_id
}

# ======================================================================
# PART 5
# QUICK LIGHTGBM TUNING GRID
#
# 10 configurations x 3 inner folds x 5 outer folds
#
# = 150 very small LightGBM fits
# ======================================================================

lightgbm_grid <- data.table::data.table(

  Config_ID = paste0(
    "LGBM_",
    sprintf(
      "%02d",
      1:10
    )
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


print(
  lightgbm_grid
)
  
# ======================================================================
# PART 6
# R -> NUMPY
# ======================================================================

to_numpy_x <- function(
  x
) {

  np$array(
    x,
    dtype = "float32"
  )
}


to_numpy_y <- function(
  y
) {

  y_numeric <- as.integer(
    as.character(
      y
    ) ==
      POSITIVE_CLASS
  )

  np$array(
    y_numeric,
    dtype = "int64"
  )
}

# ======================================================================
# PART 7
# LIGHTGBM FIT / PREDICT
# ======================================================================

fit_predict_lightgbm <- function(
  x_train,
  y_train,
  x_test,
  parameters,
  seed
) {

  x_train_py <- to_numpy_x(
    x_train
  )

  y_train_py <- to_numpy_y(
    y_train
  )

  x_test_py <- to_numpy_x(
    x_test
  )


  model <- lightgbm_py$LGBMClassifier(

    objective = "binary",

    n_estimators =
      as.integer(
        parameters$n_estimators[[1L]]
      ),

    learning_rate =
      parameters$learning_rate[[1L]],

    num_leaves =
      as.integer(
        parameters$num_leaves[[1L]]
      ),

    max_depth =
      as.integer(
        parameters$max_depth[[1L]]
      ),

    min_child_samples =
      as.integer(
        parameters$min_child_samples[[1L]]
      ),

    subsample =
      parameters$subsample[[1L]],

    subsample_freq =
      1L,

    colsample_bytree =
      parameters$colsample_bytree[[1L]],

    reg_alpha =
      parameters$reg_alpha[[1L]],

    reg_lambda =
      parameters$reg_lambda[[1L]],

    class_weight =
      "balanced",

    random_state =
      as.integer(
        seed
      ),

    n_jobs =
      -1L,

    verbosity =
      -1L
  )


  model$fit(
    x_train_py,
    y_train_py
  )


  probability_py <- model$predict_proba(
    x_test_py
  )


  probability_matrix <- reticulate::py_to_r(
    probability_py
  )


  as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )
}

# ======================================================================
# PART 8
# TabPFN-3 FOUNDATION MODEL
#
# Default pretrained model.
#
# NO gradient fine-tuning.
# NO hyperparameter search.
# NO external scaling.
# ======================================================================

# ======================================================================
# FIXED TABPFN FUNCTION FOR INSTALLED VERSION 0.1.11
#
# Installed constructor:
#
# TabPFNClassifier(
#   device = "cpu",
#   ...
#   seed = 0,
#   ...
# )
#
# IMPORTANT:
#   Version 0.1.11 uses `seed`, NOT `random_state`.
# ======================================================================

TABPFN_LABEL <- "TabPFN"

fit_predict_tabpfn <- function(
  x_train,
  y_train,
  x_test,
  seed
) {

  # ====================================================================
  # CONVERT TO NUMPY
  # ====================================================================

  x_train_py <- to_numpy_x(
    x_train
  )

  y_train_py <- to_numpy_y(
    y_train
  )

  x_test_py <- to_numpy_x(
    x_test
  )


  # ====================================================================
  # SAFE INTEGER SEED
  # ====================================================================

  seed_integer <- as.integer(
    seed %%
      .Machine$integer.max
  )


  # ====================================================================
  # SEED NUMPY
  # ====================================================================

  np$random$seed(
    seed_integer
  )


  # ====================================================================
  # CONSTRUCT TABPFN
  #
  # For version 0.1.11:
  #
  #   seed
  #
  # rather than:
  #
  #   random_state
  # ====================================================================

  model <- tabpfn_py$TabPFNClassifier(
    device = TABPFN_DEVICE,
    seed = seed_integer
  )


  # ====================================================================
  # FIT
  # ====================================================================

  model$fit(
    x_train_py,
    y_train_py
  )


  # ====================================================================
  # PREDICT PROBABILITIES
  # ====================================================================

  probability_py <- model$predict_proba(
    x_test_py
  )


  probability_matrix <- reticulate::py_to_r(
    probability_py
  )


  # ====================================================================
  # AUDIT
  # ====================================================================

  if (
    !is.matrix(
      probability_matrix
    )
  ) {

    probability_matrix <- as.matrix(
      probability_matrix
    )
  }


  if (
    ncol(
      probability_matrix
    ) !=
      2L
  ) {

    stop(
      "Expected two TabPFN probability columns, received ",
      ncol(
        probability_matrix
      ),
      "."
    )
  }


  if (
    any(
      !is.finite(
        probability_matrix
      )
    )
  ) {

    stop(
      "TabPFN returned non-finite probabilities."
    )
  }


  # ====================================================================
  # POSITIVE CLASS
  #
  # to_numpy_y() encoded:
  #
  # control = 0
  # case    = 1
  #
  # therefore column 2 = P(case)
  # ====================================================================

  probability_case <- as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )


  if (
    any(
      probability_case <
        0 |
        probability_case >
          1
    )
  ) {

    stop(
      "TabPFN probabilities fall outside [0,1]."
    )
  }


  probability_case
}

# ======================================================================
# PART 9
# PERFORMANCE METRICS
# ======================================================================

evaluate_binary_predictions <- function(
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


  auroc <- as.numeric(
    pROC::auc(
      roc_object
    )
  )


  pr_object <- PRROC::pr.curve(

    scores.class0 =
      probability[
        truth ==
          POSITIVE_CLASS
      ],

    scores.class1 =
      probability[
        truth ==
          NEGATIVE_CLASS
      ],

    curve =
      FALSE
  )


  response <- ifelse(
    probability >= threshold,
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


  data.table::data.table(

    AUROC =
      auroc,

    AUPRC =
      pr_object$auc.integral,

    BACC =
      mean(
        c(
          sensitivity,
          specificity
        )
      ),

    Sensitivity =
      sensitivity,

    Specificity =
      specificity,

    Brier =
      mean(
        (
          probability -
            as.integer(
              truth ==
                POSITIVE_CLASS
            )
        )^2
      )
  )
}

# ======================================================================
# PART 10
# SHARED OUTER FOLDS
# ======================================================================

outer_fold_id <- make_joint_stratified_folds(
  metadata = development_metadata,
  number_of_folds = OUTER_FOLDS,
  seed = SEED_FOUNDATION
)


development_metadata[
  ,
  foundation_outer_fold :=
    outer_fold_id
]


cat(
  "\nOuter fold audit:\n"
)


print(
  development_metadata[
    ,
    .(
      n = .N,
      Cases = sum(
        outcome ==
          POSITIVE_CLASS
      ),
      Controls = sum(
        outcome ==
          NEGATIVE_CLASS
      )
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
)

# ======================================================================
# PART 11
# NESTED CV
# ======================================================================

outer_prediction_rows <- list()

outer_metric_rows <- list()

inner_tuning_rows <- list()

selected_gene_rows <- list()


prediction_counter <- 0L

metric_counter <- 0L

tuning_counter <- 0L

gene_counter <- 0L


for (
  outer_fold in seq_len(
    OUTER_FOLDS
  )
) {

  cat(
    "\n\n============================================================\n",
    "OUTER FOLD ",
    outer_fold,
    " OF ",
    OUTER_FOLDS,
    "\n",
    "============================================================\n",
    sep = ""
  )


  # ====================================================================
  # OUTER TRAIN / TEST
  # ====================================================================

  outer_train_metadata <- development_metadata[
    foundation_outer_fold !=
      outer_fold
  ]


  outer_test_metadata <- development_metadata[
    foundation_outer_fold ==
      outer_fold
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


  # ====================================================================
  # INNER FOLDS
  # ====================================================================

  inner_fold_id <- make_joint_stratified_folds(
    metadata = outer_train_metadata,
    number_of_folds = INNER_FOLDS,
    seed = SEED_FOUNDATION +
      outer_fold *
        100L
  )


  outer_train_metadata[
    ,
    foundation_inner_fold :=
      inner_fold_id
  ]


  # ====================================================================
  # TUNE LIGHTGBM
  #
  # Feature selection is repeated inside EACH inner split.
  # ====================================================================

  candidate_inner_rows <- list()

  candidate_inner_counter <- 0L


  for (
    inner_fold in seq_len(
      INNER_FOLDS
    )
  ) {

    message(
      "  Inner fold ",
      inner_fold,
      "/",
      INNER_FOLDS
    )


    inner_train_metadata <- outer_train_metadata[
      foundation_inner_fold !=
        inner_fold
    ]


    inner_valid_metadata <- outer_train_metadata[
      foundation_inner_fold ==
        inner_fold
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


    # ================================================================
    # logCPM
    # ================================================================

    inner_representation <- fit_apply_logcpm(

      train_counts =
        inner_train_counts,

      new_counts =
        inner_valid_counts,

      train_metadata =
        inner_train_metadata
    )


    # ================================================================
    # TRAINING-ONLY FEATURE SELECTION
    # ================================================================

    inner_genes <- select_rf_features_final(

      train_matrix =
        inner_representation$train,

      train_metadata =
        inner_train_metadata,

      top_n =
        FOUNDATION_TOP_N,

      seed =
        SEED_FOUNDATION +
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


    # ================================================================
    # EACH LIGHTGBM CONFIG
    # ================================================================

    for (
      config_index in seq_len(
        nrow(
          lightgbm_grid
        )
      )
    ) {

      current_config <- lightgbm_grid[
        config_index
      ]


      current_probability <- fit_predict_lightgbm(

        x_train =
          inner_x_train,

        y_train =
          inner_train_metadata$outcome,

        x_test =
          inner_x_valid,

        parameters =
          current_config,

        seed =
          SEED_FOUNDATION +
          outer_fold *
            100000L +
          inner_fold *
            1000L +
          config_index
      )


      current_metric <- evaluate_binary_predictions(

        truth =
          inner_valid_metadata$outcome,

        probability =
          current_probability
      )


      candidate_inner_counter <-
        candidate_inner_counter +
        1L


      candidate_inner_rows[[candidate_inner_counter]] <- data.table::data.table(

        Outer_fold =
          outer_fold,

        Inner_fold =
          inner_fold,

        Config_ID =
          current_config$Config_ID[[1L]],

        AUROC =
          current_metric$AUROC[[1L]]
      )
    }
  }


  candidate_inner_dt <- data.table::rbindlist(
    candidate_inner_rows
  )


  candidate_summary <- candidate_inner_dt[
    ,
    .(
      Mean_inner_AUROC =
        mean(
          AUROC
        ),

      SD_inner_AUROC =
        stats::sd(
          AUROC
        )
    ),
    by =
      Config_ID
  ]


  data.table::setorder(
    candidate_summary,
    -Mean_inner_AUROC,
    SD_inner_AUROC,
    Config_ID
  )


  best_config_id <- candidate_summary$Config_ID[[1L]]


  best_lightgbm_config <- lightgbm_grid[
    Config_ID ==
      best_config_id
  ]


  cat(
    "\nSelected LightGBM config: ",
    best_config_id,
    "\n",
    "Mean inner AUROC: ",
    sprintf(
      "%.3f",
      candidate_summary$Mean_inner_AUROC[[1L]]
    ),
    "\n",
    sep = ""
  )


  tuning_counter <- tuning_counter +
    1L


  inner_tuning_rows[[tuning_counter]] <-
    data.table::copy(
      candidate_summary
    )[
      ,
      Outer_fold :=
        outer_fold
    ]


  # ====================================================================
  # OUTER FEATURE SELECTION
  #
  # Same features supplied to LightGBM AND TabPFN.
  # ====================================================================

  outer_representation <- fit_apply_logcpm(

    train_counts =
      outer_train_counts,

    new_counts =
      outer_test_counts,

    train_metadata =
      outer_train_metadata
  )


  outer_genes <- select_rf_features_final(

    train_matrix =
      outer_representation$train,

    train_metadata =
      outer_train_metadata,

    top_n =
      FOUNDATION_TOP_N,

    seed =
      SEED_FOUNDATION +
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


  selected_gene_rows[[gene_counter]] <-
    data.table::data.table(

      Outer_fold =
        outer_fold,

      Rank =
        seq_along(
          outer_genes
        ),

      Gene =
        outer_genes
    )


  # ====================================================================
  # LIGHTGBM OUTER PREDICTIONS
  # ====================================================================

  lightgbm_start <- Sys.time()


  lightgbm_probability <- fit_predict_lightgbm(

    x_train =
      outer_x_train,

    y_train =
      outer_train_metadata$outcome,

    x_test =
      outer_x_test,

    parameters =
      best_lightgbm_config,

    seed =
      SEED_FOUNDATION +
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


  # ====================================================================
  # TabPFN-3 OUTER PREDICTIONS
  # ====================================================================

  tabpfn_start <- Sys.time()


  tabpfn_probability <- fit_predict_tabpfn(

    x_train =
      outer_x_train,

    y_train =
      outer_train_metadata$outcome,

    x_test =
      outer_x_test,

    seed =
      SEED_FOUNDATION +
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


  # ====================================================================
  # SAVE PREDICTIONS
  # ====================================================================

  for (
    model_name in c(
      "LightGBM",
      "TabPFN-3"
    )
  ) {

    current_probability <- if (
      model_name ==
        "LightGBM"
    ) {
      lightgbm_probability
    } else {
      tabpfn_probability
    }


    current_seconds <- if (
      model_name ==
        "LightGBM"
    ) {
      lightgbm_seconds
    } else {
      tabpfn_seconds
    }


    prediction_counter <-
      prediction_counter +
      1L


    outer_prediction_rows[[prediction_counter]] <- data.table::data.table(

      Outer_fold =
        outer_fold,

      Model =
        model_name,

      sample_id =
        outer_test_metadata$sample_id,

      cohort_group =
        outer_test_metadata$cohort_group,

      truth =
        as.character(
          outer_test_metadata$outcome
        ),

      probability_case =
        current_probability
    )


    current_metrics <- evaluate_binary_predictions(

      truth =
        outer_test_metadata$outcome,

      probability =
        current_probability
    )


    metric_counter <-
      metric_counter +
      1L


    outer_metric_rows[[metric_counter]] <- data.table::data.table(

      Outer_fold =
        outer_fold,

      Model =
        model_name,

      AUROC =
        current_metrics$AUROC[[1L]],

      AUPRC =
        current_metrics$AUPRC[[1L]],

      BACC =
        current_metrics$BACC[[1L]],

      Brier =
        current_metrics$Brier[[1L]],

      Runtime_seconds =
        current_seconds,

      LightGBM_Config =
        ifelse(
          model_name ==
            "LightGBM",
          best_config_id,
          NA_character_
        )
    )
  }
}

# ======================================================================
# PART 12
# COMBINE RESULTS
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


data.table::fwrite(
  foundation_predictions,
  file.path(
    FOUNDATION_DIR,
    "01_nestedCV_predictions.csv"
  )
)


data.table::fwrite(
  foundation_fold_metrics,
  file.path(
    FOUNDATION_DIR,
    "02_nestedCV_fold_metrics.csv"
  )
)


data.table::fwrite(
  foundation_inner_tuning,
  file.path(
    FOUNDATION_DIR,
    "03_LightGBM_inner_tuning.csv"
  )
)


data.table::fwrite(
  foundation_selected_genes,
  file.path(
    FOUNDATION_DIR,
    "04_outer_fold_selected_genes.csv"
  )
)

# ======================================================================
# PART 13
# FOLD SUMMARY
# ======================================================================

foundation_fold_summary <- foundation_fold_metrics[
  ,
  .(

    Outer_folds =
      .N,

    Mean_AUROC =
      mean(
        AUROC
      ),

    SD_AUROC =
      stats::sd(
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
        AUPRC
      ),

    Mean_BACC =
      mean(
        BACC
      ),

    Mean_Brier =
      mean(
        Brier
      ),

    Total_runtime_seconds =
      sum(
        Runtime_seconds
      )
  ),
  by =
    Model
]


data.table::setorder(
  foundation_fold_summary,
  -Mean_AUROC
)


print(
  foundation_fold_summary
)


data.table::fwrite(
  foundation_fold_summary,
  file.path(
    FOUNDATION_DIR,
    "05_model_fold_summary.csv"
  )
)

# ======================================================================
# PART 14
# POOLED OUT-OF-FOLD PERFORMANCE
# ======================================================================

model_names <- c(
  "LightGBM",
  "TabPFN-3"
)


pooled_rows <- list()


for (
  model_index in seq_along(
    model_names
  )
) {

  current_model <- model_names[[model_index]]


  current_dt <- foundation_predictions[
    Model ==
      current_model
  ]


  stopifnot(
    data.table::uniqueN(
      current_dt$sample_id
    ) ==
      nrow(
        development_metadata
      )
  )


  current_roc <- pROC::roc(

    response = factor(
      current_dt$truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),

    predictor =
      current_dt$probability_case,

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

    scores.class0 =
      current_dt$probability_case[
        current_dt$truth ==
          POSITIVE_CLASS
      ],

    scores.class1 =
      current_dt$probability_case[
        current_dt$truth ==
          NEGATIVE_CLASS
      ],

    curve = FALSE
  )


  current_brier <- mean(
    (
      current_dt$probability_case -
        as.integer(
          current_dt$truth ==
            POSITIVE_CLASS
        )
    )^2
  )


  pooled_rows[[model_index]] <-
    data.table::data.table(

      Model =
        current_model,

      n =
        nrow(
          current_dt
        ),

      AUROC =
        as.numeric(
          pROC::auc(
            current_roc
          )
        ),

      AUROC_lower95 =
        current_ci[[1L]],

      AUROC_upper95 =
        current_ci[[3L]],

      AUPRC =
        current_pr$auc.integral,

      Brier =
        current_brier
    )
}


foundation_pooled_summary <- data.table::rbindlist(
  pooled_rows
)


print(
  foundation_pooled_summary
)


data.table::fwrite(
  foundation_pooled_summary,
  file.path(
    FOUNDATION_DIR,
    "06_pooled_OOF_performance.csv"
  )
)

# ======================================================================
# PART 15
# PAIRED MODEL COMPARISON
# ======================================================================

lightgbm_oof <- foundation_predictions[
  Model ==
    "LightGBM"
][
  order(
    sample_id
  )
]


tabpfn_oof <- foundation_predictions[
  Model ==
    "TabPFN-3"
][
  order(
    sample_id
  )
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
  pROC::auc(
    lightgbm_roc
  )
)


tabpfn_auc <- as.numeric(
  pROC::auc(
    tabpfn_roc
  )
)


foundation_comparison <- data.table::data.table(

  LightGBM_AUROC =
    lightgbm_auc,

  TabPFN3_AUROC =
    tabpfn_auc,

  Delta_AUROC_TabPFN_minus_LightGBM =
    tabpfn_auc -
    lightgbm_auc,

  Paired_DeLong_p =
    as.numeric(
      paired_delong$p.value
    )
)


print(
  foundation_comparison
)


data.table::fwrite(
  foundation_comparison,
  file.path(
    FOUNDATION_DIR,
    "07_paired_LightGBM_vs_TabPFN.csv"
  )
)

# ======================================================================
# PART 16
# PUBLICATION TABLE
# ======================================================================

publication_foundation_table <- merge(

  foundation_pooled_summary,

  foundation_fold_summary[
    ,
    .(
      Model,
      Mean_outer_AUROC =
        Mean_AUROC,
      SD_outer_AUROC =
        SD_AUROC,
      Worst_outer_AUROC =
        Worst_AUROC,
      Runtime_seconds =
        Total_runtime_seconds
    )
  ],

  by =
    "Model"
)


publication_foundation_table[
  ,
  `AUROC (95% CI)` :=
    sprintf(
      "%.3f (%.3f\u2013%.3f)",
      AUROC,
      AUROC_lower95,
      AUROC_upper95
    )
]


publication_foundation_table[
  ,
  `Mean outer AUROC (SD)` :=
    sprintf(
      "%.3f (%.3f)",
      Mean_outer_AUROC,
      SD_outer_AUROC
    )
]


publication_foundation_table[
  ,
  `Worst outer AUROC` :=
    sprintf(
      "%.3f",
      Worst_outer_AUROC
    )
]


publication_foundation_table[
  ,
  AUPRC :=
    round(
      AUPRC,
      3L
    )
]


publication_foundation_table[
  ,
  Brier :=
    round(
      Brier,
      3L
    )
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


print(
  publication_foundation_table
)


data.table::fwrite(
  publication_foundation_table,
  file.path(
    FOUNDATION_DIR,
    "Table_LightGBM_vs_TabPFN3.csv"
  )
)

# ======================================================================
# PART 17
# QUICK COMPARISON FIGURE
# ======================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)


foundation_plot_dt <- data.table::copy(
  foundation_fold_metrics
)


foundation_plot_dt[
  ,
  Model := factor(
    Model,
    levels = c(
      "LightGBM",
      "TabPFN-3"
    )
  )
]


p_foundation_compare <- ggplot2::ggplot(
  foundation_plot_dt,
  ggplot2::aes(
    x = Model,
    y = AUROC,
    group = Outer_fold
  )
) +

  ggplot2::geom_line(
    colour = "grey65",
    linewidth = 0.65
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      shape = Model
    ),
    size = 3.2
  ) +

  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey50"
  ) +

  ggplot2::coord_cartesian(
    ylim = c(
      0.45,
      1
    )
  ) +

  ggplot2::labs(
    title = "Gradient boosting versus a tabular foundation model",
    subtitle = "Paired outer-fold AUROC from pooled nested cross-validation",
    x = NULL,
    y = "Outer-fold AUROC"
  ) +

  ggplot2::theme_bw(
    base_size = 11
  ) +

  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )


print(
  p_foundation_compare
)


ggplot2::ggsave(
  filename = file.path(
    FOUNDATION_DIR,
    "Figure_LightGBM_vs_TabPFN3_outer_AUROC.pdf"
  ),
  plot = p_foundation_compare,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    FOUNDATION_DIR,
    "Figure_LightGBM_vs_TabPFN3_outer_AUROC.png"
  ),
  plot = p_foundation_compare,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ======================================================================
# PART 18A
# SELECT THE FINAL LIGHTGBM CONFIG FROM INNER-CV RESULTS ONLY
# ======================================================================

full_lgbm_tuning_summary <- foundation_inner_tuning[
  ,
  .(
    Mean_nested_inner_AUROC =
      mean(
        Mean_inner_AUROC
      ),

    SD_nested_inner_AUROC =
      stats::sd(
        Mean_inner_AUROC
      )
  ),
  by =
    Config_ID
]


data.table::setorder(
  full_lgbm_tuning_summary,
  -Mean_nested_inner_AUROC,
  SD_nested_inner_AUROC,
  Config_ID
)


FINAL_LGBM_CONFIG_ID <-
  full_lgbm_tuning_summary$Config_ID[[1L]]


FINAL_LGBM_CONFIG <- lightgbm_grid[
  Config_ID ==
    FINAL_LGBM_CONFIG_ID
]


print(
  full_lgbm_tuning_summary
)


cat(
  "\nFinal LightGBM config: ",
  FINAL_LGBM_CONFIG_ID,
  "\n",
  sep = ""
)

# ======================================================================
# PART 18B
# FULL DEVELOPMENT FIT
# ======================================================================

full_representation <- fit_apply_logcpm(

  train_counts =
    development_counts,

  new_counts =
    final_ice_counts,

  train_metadata =
    development_metadata
)


FOUNDATION_FINAL_GENES <- select_rf_features_final(

  train_matrix =
    full_representation$train,

  train_metadata =
    development_metadata,

  top_n =
    FOUNDATION_TOP_N,

  seed =
    SEED_FOUNDATION +
    999999L
)


full_x_train <- t(
  full_representation$train[
    FOUNDATION_FINAL_GENES,
    ,
    drop = FALSE
  ]
)


full_x_ice <- t(
  full_representation$new[
    FOUNDATION_FINAL_GENES,
    ,
    drop = FALSE
  ]
)


# ======================================================================
# LIGHTGBM
# ======================================================================

final_lgbm_ice_probability <- fit_predict_lightgbm(

  x_train =
    full_x_train,

  y_train =
    development_metadata$outcome,

  x_test =
    full_x_ice,

  parameters =
    FINAL_LGBM_CONFIG,

  seed =
    SEED_FOUNDATION +
    10000001L
)


# ======================================================================
# TabPFN-3
# ======================================================================

final_tabpfn_ice_probability <- fit_predict_tabpfn(

  x_train =
    full_x_train,

  y_train =
    development_metadata$outcome,

  x_test =
    full_x_ice,

  seed =
    SEED_FOUNDATION +
    10000002L
)

# ======================================================================
# PART 18C
# DIAGNOSTIC ICE COMPARISON
# ======================================================================

diagnostic_model_rows <- list()


for (
  current_model in c(
    "LightGBM",
    "TabPFN-3"
  )
) {

  current_probability <- if (
    current_model ==
      "LightGBM"
  ) {
    final_lgbm_ice_probability
  } else {
    final_tabpfn_ice_probability
  }


  current_roc <- pROC::roc(
    response = final_ice_metadata$outcome,
    predictor = current_probability,
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
      method = "delong"
    )
  )


  current_metrics <- evaluate_binary_predictions(
    truth = final_ice_metadata$outcome,
    probability = current_probability
  )


  diagnostic_model_rows[[current_model]] <-
    data.table::data.table(

      Model =
        current_model,

      AUROC =
        as.numeric(
          pROC::auc(
            current_roc
          )
        ),

      AUROC_lower95 =
        current_ci[[1L]],

      AUROC_upper95 =
        current_ci[[3L]],

      AUPRC =
        current_metrics$AUPRC[[1L]],

      BACC =
        current_metrics$BACC[[1L]],

      Brier =
        current_metrics$Brier[[1L]]
    )
}


foundation_ice_summary <- data.table::rbindlist(
  diagnostic_model_rows
)


print(
  foundation_ice_summary
)


data.table::fwrite(
  foundation_ice_summary,
  file.path(
    FOUNDATION_DIR,
    "08_diagnostic_ICE_LightGBM_vs_TabPFN3.csv"
  )
)
