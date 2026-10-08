# ======================================================================
# PART 19
# FINAL DEVELOPMENT FIT:
# LIGHTGBM vs TABPFN
#
# TRAIN:
#   372 development samples
#
# TEST:
#   93 diagnostic final-ICE samples
#
# REPRESENTATION:
#   log2(CPM + 1)
#
# FEATURES:
#   SAME 100 RF-selected genes for both models
#
# LIGHTGBM:
#   hyperparameters selected from nested-CV inner results only
#
# TABPFN:
#   no task-specific tuning
#
# IMPORTANT:
#   ICE is NOT used for model selection or feature selection.
# ======================================================================


# ======================================================================
# 19.1 SETTINGS
# ======================================================================

ICE_COMPARISON_DIR <- file.path(
  FOUNDATION_DIR,
  "diagnostic_FINAL_ICE"
)

dir.create(
  ICE_COMPARISON_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


ICE_TOP_N <- 100L


FOUNDATION_LABEL <- if (
  exists(
    "TABPFN_LABEL"
  )
) {
  TABPFN_LABEL
} else {
  "TabPFN"
}


ICE_MODEL_ORDER <- c(
  "LightGBM",
  FOUNDATION_LABEL
)


# ======================================================================
# 19.2 SELECT FINAL LIGHTGBM CONFIGURATION
#
# Uses ONLY inner-CV results accumulated during nested CV.
# ICE plays absolutely no role.
# ======================================================================

final_lgbm_tuning_summary <- foundation_inner_tuning[
  ,
  .(

    Outer_folds =
      .N,

    Mean_inner_AUROC =
      mean(
        Mean_inner_AUROC,
        na.rm = TRUE
      ),

    SD_inner_AUROC =
      stats::sd(
        Mean_inner_AUROC,
        na.rm = TRUE
      ),

    Worst_outer_inner_AUROC =
      min(
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


FINAL_LGBM_CONFIG_ID <-
  final_lgbm_tuning_summary$Config_ID[[1L]]


FINAL_LGBM_CONFIG <- lightgbm_grid[
  Config_ID ==
    FINAL_LGBM_CONFIG_ID
]


cat(
  "\n============================================================\n",
  "FINAL LIGHTGBM CONFIGURATION\n",
  "============================================================\n",
  "Config: ",
  FINAL_LGBM_CONFIG_ID,
  "\n",
  "Mean inner AUROC: ",
  sprintf(
    "%.3f",
    final_lgbm_tuning_summary$Mean_inner_AUROC[[1L]]
  ),
  "\n",
  "SD inner AUROC: ",
  sprintf(
    "%.3f",
    final_lgbm_tuning_summary$SD_inner_AUROC[[1L]]
  ),
  "\n",
  "============================================================\n",
  sep = ""
)


print(
  FINAL_LGBM_CONFIG
)


data.table::fwrite(
  final_lgbm_tuning_summary,
  file.path(
    ICE_COMPARISON_DIR,
    "01_final_LightGBM_config_selection.csv"
  )
)

# ======================================================================
# PART 20
# TRAINING-ONLY REPRESENTATION + FEATURE SELECTION
# ======================================================================


# ======================================================================
# 20.1 logCPM
#
# Training information comes from development data only.
# ======================================================================

ice_representation <- fit_apply_logcpm(

  train_counts =
    development_counts,

  new_counts =
    final_ice_counts,

  train_metadata =
    development_metadata
)


# ======================================================================
# 20.2 FEATURE SELECTION
#
# Selection is performed ONCE using all 372 development observations.
# Exactly the same 100 genes are supplied to both models.
# ======================================================================

ICE_SELECTED_GENES <- select_rf_features_final(

  train_matrix =
    ice_representation$train,

  train_metadata =
    development_metadata,

  top_n =
    ICE_TOP_N,

  seed =
    SEED_FOUNDATION +
      99000001L
)


cat(
  "\nSelected genes: ",
  length(
    ICE_SELECTED_GENES
  ),
  "\n",
  sep = ""
)


print(
  ICE_SELECTED_GENES
)


data.table::fwrite(
  data.table::data.table(
    Rank =
      seq_along(
        ICE_SELECTED_GENES
      ),
    Gene =
      ICE_SELECTED_GENES
  ),
  file.path(
    ICE_COMPARISON_DIR,
    "02_final_100_selected_genes.csv"
  )
)


# ======================================================================
# 20.3 MODEL MATRICES
# ======================================================================

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
  nrow(
    ice_x_train
  ) ==
    nrow(
      development_metadata
    )
)


stopifnot(
  nrow(
    ice_x_test
  ) ==
    nrow(
      final_ice_metadata
    )
)


stopifnot(
  identical(
    rownames(
      ice_x_test
    ),
    final_ice_metadata$sample_id
  )
)

# ======================================================================
# PART 21
# FIT BOTH MODELS ON ALL DEVELOPMENT DATA
# ======================================================================


# ======================================================================
# 21.1 LIGHTGBM
# ======================================================================

cat(
  "\nFitting final LightGBM...\n"
)


lgbm_start_time <- Sys.time()


final_ice_lgbm_probability <- fit_predict_lightgbm(

  x_train =
    ice_x_train,

  y_train =
    development_metadata$outcome,

  x_test =
    ice_x_test,

  parameters =
    FINAL_LGBM_CONFIG,

  seed =
    SEED_FOUNDATION +
      99000002L
)


lgbm_runtime_seconds <- as.numeric(
  difftime(
    Sys.time(),
    lgbm_start_time,
    units = "secs"
  )
)


# ======================================================================
# 21.2 TABPFN
# ======================================================================

cat(
  "\nFitting final ",
  FOUNDATION_LABEL,
  "...\n",
  sep = ""
)


tabpfn_start_time <- Sys.time()


final_ice_tabpfn_probability <- fit_predict_tabpfn(

  x_train =
    ice_x_train,

  y_train =
    development_metadata$outcome,

  x_test =
    ice_x_test,

  seed =
    SEED_FOUNDATION +
      99000003L
)


tabpfn_runtime_seconds <- as.numeric(
  difftime(
    Sys.time(),
    tabpfn_start_time,
    units = "secs"
  )
)


# ======================================================================
# 21.3 HARD PROBABILITY AUDIT
# ======================================================================

stopifnot(
  length(
    final_ice_lgbm_probability
  ) ==
    nrow(
      final_ice_metadata
    )
)


stopifnot(
  length(
    final_ice_tabpfn_probability
  ) ==
    nrow(
      final_ice_metadata
    )
)


stopifnot(
  all(
    is.finite(
      final_ice_lgbm_probability
    )
  )
)


stopifnot(
  all(
    is.finite(
      final_ice_tabpfn_probability
    )
  )
)


cat(
  "\nLightGBM probability range: ",
  sprintf(
    "%.3f",
    min(
      final_ice_lgbm_probability
    )
  ),
  " - ",
  sprintf(
    "%.3f",
    max(
      final_ice_lgbm_probability
    )
  ),
  "\n",
  FOUNDATION_LABEL,
  " probability range: ",
  sprintf(
    "%.3f",
    min(
      final_ice_tabpfn_probability
    )
  ),
  " - ",
  sprintf(
    "%.3f",
    max(
      final_ice_tabpfn_probability
    )
  ),
  "\n",
  sep = ""
)
  # ======================================================================
# PART 22
# PAIRED DIAGNOSTIC ICE PREDICTIONS
# ======================================================================

ice_model_predictions <- data.table::rbindlist(
  list(

    data.table::data.table(

      Model =
        "LightGBM",

      sample_id =
        final_ice_metadata$sample_id,

      cohort_group =
        final_ice_metadata$cohort_group,

      truth =
        as.character(
          final_ice_metadata$outcome
        ),

      probability_case =
        final_ice_lgbm_probability
    ),

    data.table::data.table(

      Model =
        FOUNDATION_LABEL,

      sample_id =
        final_ice_metadata$sample_id,

      cohort_group =
        final_ice_metadata$cohort_group,

      truth =
        as.character(
          final_ice_metadata$outcome
        ),

      probability_case =
        final_ice_tabpfn_probability
    )
  )
)


ice_model_predictions[
  ,
  Model := factor(
    Model,
    levels = ICE_MODEL_ORDER
  )
]


# ======================================================================
# AUDIT
# ======================================================================

ice_prediction_audit <- ice_model_predictions[
  ,
  .(
    Rows = .N,
    Unique_samples = data.table::uniqueN(
      sample_id
    ),
    Cases = sum(
      truth ==
        POSITIVE_CLASS
    ),
    Controls = sum(
      truth ==
        NEGATIVE_CLASS
    )
  ),
  by = Model
]


print(
  ice_prediction_audit
)


data.table::fwrite(
  ice_model_predictions,
  file.path(
    ICE_COMPARISON_DIR,
    "03_FINAL_ICE_predictions_LightGBM_vs_TabPFN.csv"
  )
)

# ======================================================================
# PART 23
# MODEL PERFORMANCE ON DIAGNOSTIC FINAL ICE
# ======================================================================

ICE_THRESHOLD <- 0.50


ice_metric_rows <- list()


for (
  model_index in seq_along(
    ICE_MODEL_ORDER
  )
) {

  current_model <- ICE_MODEL_ORDER[[model_index]]


  current_dt <- ice_model_predictions[
    Model ==
      current_model
  ]


  # ====================================================================
  # ROC
  # ====================================================================

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


  current_auc <- as.numeric(
    pROC::auc(
      current_roc
    )
  )


  current_ci <- as.numeric(
    pROC::ci.auc(
      current_roc,
      method = "delong",
      conf.level = 0.95
    )
  )


  # ====================================================================
  # PR
  # ====================================================================

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


  # ====================================================================
  # THRESHOLD METRICS
  # ====================================================================

  current_response <- ifelse(
    current_dt$probability_case >=
      ICE_THRESHOLD,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  )


  current_sensitivity <- mean(
    current_response[
      current_dt$truth ==
        POSITIVE_CLASS
    ] ==
      POSITIVE_CLASS
  )


  current_specificity <- mean(
    current_response[
      current_dt$truth ==
        NEGATIVE_CLASS
    ] ==
      NEGATIVE_CLASS
  )


  current_bacc <- mean(
    c(
      current_sensitivity,
      current_specificity
    )
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


  ice_metric_rows[[model_index]] <- data.table::data.table(

    Model =
      current_model,

    n =
      nrow(
        current_dt
      ),

    AUROC =
      current_auc,

    AUROC_lower95 =
      current_ci[[1L]],

    AUROC_upper95 =
      current_ci[[3L]],

    AUPRC =
      current_pr$auc.integral,

    Prevalence =
      mean(
        current_dt$truth ==
          POSITIVE_CLASS
      ),

    BACC =
      current_bacc,

    Sensitivity =
      current_sensitivity,

    Specificity =
      current_specificity,

    Brier =
      current_brier,

    Threshold =
      ICE_THRESHOLD,

    Runtime_seconds =
      if (
        current_model ==
          "LightGBM"
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


print(
  ice_model_summary
)


data.table::fwrite(
  ice_model_summary,
  file.path(
    ICE_COMPARISON_DIR,
    "04_FINAL_ICE_performance_LightGBM_vs_TabPFN.csv"
  )
)

# ======================================================================
# PART 24
# PAIRED DELONG TEST
# ======================================================================

paired_ice_dt <- data.table::dcast(
  ice_model_predictions,
  sample_id +
    cohort_group +
    truth ~ Model,
  value.var = "probability_case"
)


# Find the generated TabPFN column safely.
foundation_probability_column <- setdiff(
  colnames(
    paired_ice_dt
  ),
  c(
    "sample_id",
    "cohort_group",
    "truth",
    "LightGBM"
  )
)


if (
  length(
    foundation_probability_column
  ) !=
    1L
) {

  stop(
    "Could not uniquely identify the TabPFN probability column."
  )
}


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

