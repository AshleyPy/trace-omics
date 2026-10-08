# ======================================================================
# PART 20A
# DEPLOYMENT PREDICTION FUNCTION
# ======================================================================

predict_tb_final_rescue <- function(
  bundle,
  new_counts
) {

  if (
    is.null(
      rownames(
        new_counts
      )
    )
  ) {

    stop(
      "new_counts requires gene row names."
    )
  }

  if (
    is.null(
      colnames(
        new_counts
      )
    )
  ) {

    stop(
      "new_counts requires sample column names."
    )
  }

  missing_genes <- setdiff(
    bundle$raw_gene_order,
    rownames(
      new_counts
    )
  )

  if (
    length(
      missing_genes
    ) >
      0L
  ) {

    stop(
      "New counts are missing ",
      length(
        missing_genes
      ),
      " genes required by the model."
    )
  }

  new_counts <- new_counts[
    bundle$raw_gene_order,
    ,
    drop = FALSE
  ]

  library_size <- colSums(
    new_counts
  )

  if (
    any(
      library_size <=
        0
    )
  ) {

    stop(
      "Invalid RNA-seq library size detected."
    )
  }

  cpm <- sweep(
    new_counts,
    MARGIN = 2L,
    STATS = library_size,
    FUN = "/"
  ) *
    1e6

  logcpm <- log2(
    cpm +
      1
  )

  x_new <- t(
    logcpm[
      bundle$selected_genes,
      ,
      drop = FALSE
    ]
  )

  probability_matrix <- predict(
    bundle$model,
    data = as.data.frame(
      x_new,
      check.names = FALSE
    )
  )$predictions

  probability_case <- as.numeric(
    probability_matrix[
      ,
      bundle$positive_class
    ]
  )

  predicted_class <- ifelse(
    probability_case >=
      bundle$threshold,
    bundle$positive_class,
    bundle$negative_class
  )

  data.table::data.table(
    sample_id = colnames(
      new_counts
    ),
    probability_case = probability_case,
    predicted_class = predicted_class
  )
}


# ======================================================================
# PART 20B
# CREATE BUNDLE
# ======================================================================

FINAL_BUNDLE <- list(

  bundle_version =
    "TB_RESCUE_FINAL_2.1",

  created =
    as.character(
      Sys.time()
    ),

  disease =
    "TB",

  input_scale =
    "raw integer gene x sample counts",

  representation =
    "log2(CPM + 1)",

  feature_selection =
    "RF permutation importance",

  selected_genes =
    FINAL_GENES,

  top_n =
    LOCKED_TOP_N,

  model_name =
    "Random Forest",

  class_weighting =
    LOCKED_WEIGHTING,

  threshold =
    LOCKED_THRESHOLD,

  negative_class =
    NEGATIVE_CLASS,

  positive_class =
    POSITIVE_CLASS,

  raw_gene_order =
    rownames(
      development_counts
    ),

  development_sample_ids =
    development_metadata$sample_id,

  strict_loco_selection =
    locked_configuration,

  model =
    FINAL_RF,

  predict_function =
    predict_tb_final_rescue
)


saveRDS(
  FINAL_BUNDLE,
  file.path(
    MODEL_DIR,
    "TB_FINAL_RESCUE_MODEL_BUNDLE.rds"
  ),
  compress = "xz"
)

# ======================================================================
# PART 21
# DEPLOYMENT SELF-TEST
# ======================================================================

bundle_self_test <- FINAL_BUNDLE$predict_function(
  bundle = FINAL_BUNDLE,
  new_counts = final_ice_counts
)


bundle_self_test <- merge(

  bundle_self_test,

  diagnostic_ice_predictions[
    ,
    .(
      sample_id,
      original_probability =
        probability_case
    )
  ],

  by = "sample_id"
)


bundle_self_test[
  ,
  Absolute_difference :=
    abs(
      probability_case -
        original_probability
    )
]


max_bundle_difference <- max(
  bundle_self_test$Absolute_difference
)


cat(
  "\nMaximum frozen-bundle probability difference: ",
  format(
    max_bundle_difference,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


if (
  max_bundle_difference >
    1e-8
) {

  stop(
    "Frozen deployment bundle does not reproduce saved predictions."
  )
}


message(
  "PASS: frozen deployment bundle reproduces saved predictions."
)

# ======================================================================
# PART 22A
# PUBLICATION MANIFEST
# ======================================================================

publication_manifest <- data.table::data.table(

  Item = c(
    "Table 1",
    "Table 2",
    "Table 3",
    "Table 4",
    "Table 5",
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5",
    "Figure 6",
    "Figure 7",
    "Figure S1",
    "Figure S2",
    "Final model bundle"
  ),

  Description = c(
    "Development-only 25/50/100-gene x class-weighting optimisation",
    "Locked classifier strict leave-one-cohort-out performance",
    "Diagnostic held-out test performance",
    "Calibration intercept, slope and Brier score",
    "Repeated permutation gene importance with 95% intervals",
    "Development-only optimisation",
    "Empirical stepwise strict LOCO ROC with pointwise 95% bootstrap bands",
    "Strict LOCO precision-recall curves",
    "Strict LOCO calibration",
    "Combined ROC, PR and calibration validation figure",
    "Post-hoc Delta-AUROC predictive-gene importance",
    "Accumulated local effect profiles",
    "Diagnostic held-out ROC, PR and calibration",
    "DALEX global feature importance",
    "Frozen logCPM + Random Forest deployment bundle"
  )
)


data.table::fwrite(
  publication_manifest,
  file.path(
    FINAL_DIR,
    "PUBLICATION_MANIFEST.csv"
  )
)


# ======================================================================
# PART 22B
# FINAL RESULT AUDIT
# ======================================================================

cat(
  "\n\n============================================================\n",
  "TB FINAL PUBLICATION MODEL COMPLETE\n",
  "============================================================\n",
  "Representation: log2(CPM + 1)\n",
  "Feature selection: RF permutation importance\n",
  "Model: Random Forest\n",
  "Genes: ",
  LOCKED_TOP_N,
  "\n",
  "Class weighting: ",
  LOCKED_WEIGHTING,
  "\n",
  "Threshold: ",
  sprintf(
    "%.2f",
    LOCKED_THRESHOLD
  ),
  "\n\n",
  "Mean strict LOCO AUROC: ",
  sprintf(
    "%.3f",
    locked_configuration$Mean_LOCO_AUROC[[1L]]
  ),
  "\n",
  "Worst strict LOCO AUROC: ",
  sprintf(
    "%.3f",
    locked_configuration$Worst_cohort_AUROC[[1L]]
  ),
  "\n\n",
  "Diagnostic ICE AUROC: ",
  sprintf(
    "%.3f",
    ice_auc
  ),
  " (",
  sprintf(
    "%.3f",
    ice_auc_ci[[1L]]
  ),
  "\u2013",
  sprintf(
    "%.3f",
    ice_auc_ci[[3L]]
  ),
  ")\n",
  "Diagnostic ICE AUPRC: ",
  sprintf(
    "%.3f",
    ice_pr_object$auc.integral
  ),
  "\n",
  "Diagnostic ICE BACC: ",
  sprintf(
    "%.3f",
    diagnostic_ice_publication_result$BACC[[1L]]
  ),
  "\n",
  "Diagnostic ICE sensitivity: ",
  sprintf(
    "%.3f",
    diagnostic_ice_publication_result$Sensitivity[[1L]]
  ),
  "\n",
  "Diagnostic ICE specificity: ",
  sprintf(
    "%.3f",
    diagnostic_ice_publication_result$Specificity[[1L]]
  ),
  "\n",
  "Bundle reproduction error: ",
  format(
    max_bundle_difference,
    scientific = TRUE
  ),
  "\n\n",
  "ROC display: EMPIRICAL STEPWISE\n",
  "AUROC CI: DeLong 95%\n",
  "ROC ribbon: pointwise stratified-bootstrap 95%\n",
  "\nPublication directory:\n",
  FINAL_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)
