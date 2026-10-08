# ============================================================================
# 21. SELECTED-GENE STABILITY ACROSS CV FOLDS
# ============================================================================

# Check required columns first
required_columns <- c(
  "Model_ID",
  "Gene"
)

missing_columns <- setdiff(
  required_columns,
  names(cv_selected_features)
)

if (length(missing_columns) > 0L) {
  stop(
    "cv_selected_features is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ---------------------------------------------------------------------------
# 21.1 Calculate selection frequency
# ---------------------------------------------------------------------------

feature_stability <- cv_selected_features[
  ,
  .(
    Selection_count = .N,
    Selection_frequency = .N / (CV_FOLDS * CV_REPEATS)
  ),
  by = .(
    Model_ID,
    Gene
  )
]


# ---------------------------------------------------------------------------
# 21.2 Add readable model labels
# ---------------------------------------------------------------------------

model_labels <- c(
  Model_A_TB_naive = "Model A: TB-naive",
  Model_B_ALL5 = "Model B: ALL5"
)

feature_stability[
  ,
  Model := unname(
    model_labels[Model_ID]
  )
]


# Check that every model ID was successfully mapped
if (anyNA(feature_stability$Model)) {

  unmapped_models <- unique(
    feature_stability[
      is.na(Model),
      Model_ID
    ]
  )

  stop(
    "The following Model_ID values have no model label: ",
    paste(unmapped_models, collapse = ", ")
  )
}


# ---------------------------------------------------------------------------
# 21.3 Reorder and rank
# ---------------------------------------------------------------------------

data.table::setcolorder(
  feature_stability,
  c(
    "Model_ID",
    "Model",
    "Gene",
    "Selection_count",
    "Selection_frequency"
  )
)

data.table::setorder(
  feature_stability,
  Model_ID,
  -Selection_frequency,
  Gene
)

feature_stability[
  ,
  Stability_rank := seq_len(.N),
  by = Model_ID
]


# ---------------------------------------------------------------------------
# 21.4 Audit
# ---------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "SELECTED-GENE STABILITY\n",
  "============================================================\n",
  sep = ""
)

print(
  feature_stability[
    Selection_frequency >= 0.50
  ]
)

# ============================================================================
# 22. MANIFEST + MASTER RESULT
# ============================================================================

manifest <- data.table(
  Item = c(
    "Pipeline",
    "Feature universe",
    "Development CV",
    "Representation",
    "FS method",
    "Learner",
    "Threshold source",
    "Strict target-batch rule",
    "Output root"
  ),
  Value = c(
    PIPELINE_LABEL,
    "Fixed 4,450 genes shared across TB, RA, UC, CD and T1D",
    paste0(CV_FOLDS, "-fold x ", CV_REPEATS, " repeats; cohort x outcome stratified"),
    "Training-fitted VST followed by training-fitted ComBat-Ref",
    paste0("Boruta top ", TOP_N, " within each CV training fold"),
    paste0("Fixed LightGBM: iterations=", LGB_NUM_ITERATIONS,
           "; min_data_in_leaf=", LGB_MIN_DATA_IN_LEAF,
           "; bagging_fraction=", LGB_BAGGING_FRACTION),
    "Youden threshold from averaged development OOF predictions only",
    "No target-batch parameter estimation; unseen target cohorts cause final evaluation to be skipped",
    candidate_root
  )
)

fwrite(manifest, file.path(candidate_root, "CANDIDATE_EXPERIMENT_MANIFEST.csv"))

master <- list(
  settings = list(
    seed = SEED,
    pipeline = PIPELINE_LABEL,
    top_n = TOP_N,
    cv_folds = CV_FOLDS,
    cv_repeats = CV_REPEATS,
    strict_mode = STRICT_MODE,
    lightgbm = list(
      num_iterations = LGB_NUM_ITERATIONS,
      min_data_in_leaf = LGB_MIN_DATA_IN_LEAF,
      bagging_fraction = LGB_BAGGING_FRACTION,
      bagging_freq = LGB_BAGGING_FREQ
    )
  ),
  input_audit = input_audit,
  final_ice_audit = final_ice_audit,
  cv_fold_metrics = cv_fold_metrics,
  cv_summary = cv_summary,
  frozen_thresholds = threshold_table,
  final_status = final_status,
  final_results = final_results,
  feature_stability = feature_stability
)

saveRDS(master, file.path(candidate_root, "MASTER_candidate_VST_ComBatRef_Boruta_LightGBM.rds"))

cat(
  "\n============================================================\n",
  "CANDIDATE PIPELINE ANALYSIS COMPLETE\n",
  "============================================================\n",
  "Output directory:\n", candidate_root, "\n\n",
  "Development CV summary:\n",
  sep = ""
)
print(cv_summary)

cat("\nFrozen thresholds:\n")
print(threshold_table)

cat("\nFinal evaluation status:\n")
print(final_status)

if (final_results[["Model_A_TB_naive"]]$status != "COMPLETE") {
  cat(
    "\nIMPORTANT:\n",
    "Model A final-ice candidate evaluation was intentionally not forced.\n",
    "The TB-naive development population lacks one or more TB final-ice cohorts,\n",
    "so cohort-level ComBat-Ref has no training estimates for those batches.\n",
    "Estimating them from final-ice data would make the analysis transductive.\n",
    "Keep the existing VST + fSVA pipeline for strict TB-naive / LODO transfer,\n",
    "or report a separately labelled transductive ComBat-Ref analysis if you\n",
    "explicitly decide that unsupervised target adaptation is part of the design.\n",
    sep = ""
  )
}


