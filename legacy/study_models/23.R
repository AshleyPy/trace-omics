# ======================================================================
# PART 23
# DEVELOPMENT-vs-TEST GENERALISABILITY TABLE
#
# "TRAIN" PERFORMANCE IS NOT APPARENT TRAINING PERFORMANCE.
#
# Development AUROC:
#   pooled strict leave-one-cohort-out predictions
#
# Test AUROC:
#   diagnostic held-out ICE predictions
#
# Both:
#   empirical AUROC
#   DeLong 95% CI
# ======================================================================


# ======================================================================
# 23.1 AUDIT DEVELOPMENT LOCO PREDICTIONS
# ======================================================================

development_generalisation_dt <- data.table::copy(
  locked_loco_predictions
)


cat(
  "\n============================================================\n",
  "DEVELOPMENT LOCO PREDICTION AUDIT\n",
  "============================================================\n",
  "Rows: ",
  nrow(
    development_generalisation_dt
  ),
  "\n",
  "Unique samples: ",
  data.table::uniqueN(
    development_generalisation_dt$sample_id
  ),
  "\n",
  "Expected development samples: ",
  nrow(
    development_metadata
  ),
  "\n",
  sep = ""
)


if (
  data.table::uniqueN(
    development_generalisation_dt$sample_id
  ) !=
    nrow(
      development_metadata
    )
) {

  stop(
    "LOCO predictions do not contain exactly one prediction per development sample."
  )
}


if (
  anyDuplicated(
    development_generalisation_dt$sample_id
  ) >
    0L
) {

  stop(
    "Duplicate development sample predictions detected."
  )
}


# ======================================================================
# 23.2 DEVELOPMENT STRICT-LOCO ROC
# ======================================================================

development_loco_roc <- pROC::roc(

  response = factor(
    development_generalisation_dt$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),

  predictor =
    development_generalisation_dt$probability_case,

  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),

  direction = "<",

  quiet = TRUE
)


development_loco_auc <- as.numeric(
  pROC::auc(
    development_loco_roc
  )
)


development_loco_ci <- as.numeric(
  pROC::ci.auc(
    development_loco_roc,
    method = "delong",
    conf.level = 0.95
  )
)


# ======================================================================
# 23.3 TEST ROC
# ======================================================================

test_roc <- pROC::roc(

  response = factor(
    diagnostic_ice_predictions$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),

  predictor =
    diagnostic_ice_predictions$probability_case,

  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),

  direction = "<",

  quiet = TRUE
)


test_auc <- as.numeric(
  pROC::auc(
    test_roc
  )
)


test_ci <- as.numeric(
  pROC::ci.auc(
    test_roc,
    method = "delong",
    conf.level = 0.95
  )
)


# ======================================================================
# 23.4 GENERALISATION GAP
#
# Positive:
#   test AUROC exceeds development LOCO AUROC
#
# Negative:
#   discrimination fell in the test population
# ======================================================================

delta_test_minus_development <- (
  test_auc -
    development_loco_auc
)


absolute_generalisation_gap <- abs(
  delta_test_minus_development
)


# ======================================================================
# 23.5 ONE-ROW PUBLICATION TABLE
# ======================================================================

generalisability_table <- data.table::data.table(

  Model =
    paste0(
      "logCPM + RF | ",
      LOCKED_TOP_N,
      " genes",
      ifelse(
        LOCKED_WEIGHTING,
        " | weighted",
        " | unweighted"
      )
    ),

  Development_n =
    nrow(
      development_generalisation_dt
    ),

  `Development strict-LOCO AUROC (95% CI)` =
    sprintf(
      "%.3f (%.3f\u2013%.3f)",
      development_loco_auc,
      development_loco_ci[[1L]],
      development_loco_ci[[3L]]
    ),

  Test_n =
    nrow(
      diagnostic_ice_predictions
    ),

  `Test AUROC (95% CI)` =
    sprintf(
      "%.3f (%.3f\u2013%.3f)",
      test_auc,
      test_ci[[1L]],
      test_ci[[3L]]
    ),

  `Delta AUROC (Test - Development)` =
    sprintf(
      "%+.3f",
      delta_test_minus_development
    ),

  `Absolute generalisation gap` =
    sprintf(
      "%.3f",
      absolute_generalisation_gap
    )
)


print(
  generalisability_table
)


# ======================================================================
# 23.6 SAVE CSV
# ======================================================================

data.table::fwrite(
  generalisability_table,
  file.path(
    TABLE_DIR,
    "Table_6_development_vs_test_generalisability.csv"
  )
)
