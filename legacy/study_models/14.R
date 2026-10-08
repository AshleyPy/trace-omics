# ======================================================================
# PART 12
# DIAGNOSTIC ICE PERFORMANCE TABLE
# ======================================================================

ice_threshold_metrics <- calculate_threshold_metrics_publication(
  truth = diagnostic_ice_predictions$truth,
  probability = diagnostic_ice_predictions$probability_case,
  threshold = LOCKED_THRESHOLD
)


diagnostic_ice_publication_result <- data.table::data.table(

  Representation = "log2(CPM + 1)",

  Feature_selection =
    "RF permutation importance",

  TOP_N =
    LOCKED_TOP_N,

  Model =
    "Random Forest",

  Class_weighting =
    LOCKED_WEIGHTING,

  Threshold =
    LOCKED_THRESHOLD,

  n =
    nrow(
      diagnostic_ice_predictions
    ),

  Cases =
    sum(
      diagnostic_ice_predictions$truth ==
        POSITIVE_CLASS
    ),

  Controls =
    sum(
      diagnostic_ice_predictions$truth ==
        NEGATIVE_CLASS
    ),

  AUROC =
    ice_auc,

  AUROC_lower95 =
    ice_auc_ci[[1L]],

  AUROC_upper95 =
    ice_auc_ci[[3L]],

  AUPRC =
    ice_pr_object$auc.integral,

  Prevalence =
    ice_prevalence,

  BACC =
    ice_threshold_metrics$BACC[[1L]],

  Sensitivity =
    ice_threshold_metrics$Sensitivity[[1L]],

  Specificity =
    ice_threshold_metrics$Specificity[[1L]],

  Brier =
    ice_threshold_metrics$Brier[[1L]]
)


print(
  diagnostic_ice_publication_result
)


data.table::fwrite(
  diagnostic_ice_publication_result,
  file.path(
    TABLE_DIR,
    "Table_3_diagnostic_ICE_performance.csv"
  )
)
