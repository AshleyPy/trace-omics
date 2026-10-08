ice_tabpfn_roc <- pROC::roc(
  response = paired_ice_truth,
  predictor = paired_ice_dt[[foundation_probability_column]],
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


ice_paired_delong <- pROC::roc.test(
  ice_lgbm_roc,
  ice_tabpfn_roc,
  method = "delong",
  paired = TRUE
)


ice_lgbm_auc <- as.numeric(
  pROC::auc(
    ice_lgbm_roc
  )
)


ice_tabpfn_auc <- as.numeric(
  pROC::auc(
    ice_tabpfn_roc
  )
)


ice_delta_auc <- (
  ice_tabpfn_auc -
    ice_lgbm_auc
)


ice_delong_p <- as.numeric(
  ice_paired_delong$p.value
)


ice_paired_comparison <- data.table::data.table(

  LightGBM_AUROC =
    ice_lgbm_auc,

  TabPFN_AUROC =
    ice_tabpfn_auc,

  Delta_AUROC_TabPFN_minus_LightGBM =
    ice_delta_auc,

  Paired_DeLong_p =
    ice_delong_p
)


print(
  ice_paired_comparison
)


data.table::fwrite(
  ice_paired_comparison,
  file.path(
    ICE_COMPARISON_DIR,
    "05_FINAL_ICE_paired_DeLong_comparison.csv"
  )
)

# ======================================================================
# PART 25
# PUBLICATION TABLE
# ======================================================================

ice_publication_table <- data.table::copy(
  ice_model_summary
)


ice_publication_table[
  ,
  `AUROC (95% CI)` := sprintf(
    "%.3f (%.3f\u2013%.3f)",
    AUROC,
    AUROC_lower95,
    AUROC_upper95
  )
]


ice_publication_table[
  ,
  `:=`(
    AUPRC =
      round(
        AUPRC,
        3L
      ),
    BACC =
      round(
        BACC,
        3L
      ),
    Sensitivity =
      round(
        Sensitivity,
        3L
      ),
    Specificity =
      round(
        Specificity,
        3L
      ),
    Brier =
      round(
        Brier,
        3L
      ),
    Runtime_seconds =
      round(
        Runtime_seconds,
        2L
      )
  )
]


ice_publication_table <- ice_publication_table[
  ,
  .(
    Model,
    n,
    `AUROC (95% CI)`,
    AUPRC,
    BACC,
    Sensitivity,
    Specificity,
    Brier,
    Runtime_seconds
  )
]


print(
  ice_publication_table
)


data.table::fwrite(
  ice_publication_table,
  file.path(
    ICE_COMPARISON_DIR,
    "Table_FINAL_ICE_LightGBM_vs_TabPFN.csv"
  )
)
