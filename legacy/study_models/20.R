# ======================================================================
# PART 18A
# TABLE 1
# SMALL MODEL OPTIMISATION
# ======================================================================

publication_optimisation_table <- optimisation_summary[
  ,
  .(

    Rank,

    Genes =
      TOP_N,

    `Class weighting` =
      ifelse(
        Class_weighting,
        "Yes",
        "No"
      ),

    `Mean LOCO AUROC` =
      round(
        Mean_LOCO_AUROC,
        3L
      ),

    `Worst cohort AUROC` =
      round(
        Worst_cohort_AUROC,
        3L
      ),

    `SD AUROC` =
      round(
        SD_LOCO_AUROC,
        3L
      ),

    `Mean LOCO AUPRC` =
      round(
        Mean_LOCO_AUPRC,
        3L
      ),

    `Mean BACC` =
      round(
        Mean_LOCO_BACC,
        3L
      ),

    `Mean Brier` =
      round(
        Mean_Brier,
        3L
      )
  )
]


print(
  publication_optimisation_table
)


data.table::fwrite(
  publication_optimisation_table,
  file.path(
    TABLE_DIR,
    "Table_1_model_optimisation.csv"
  )
)


# ======================================================================
# PART 18B
# TABLE 2
# LOCKED MODEL LOCO PERFORMANCE
# ======================================================================

locked_loco_metrics <- optimisation_metrics[
  Config_ID ==
    LOCKED_CONFIG_ID
]


publication_loco_table <- locked_loco_metrics[
  ,
  .(

    Cohort =
      Target_cohort,

    n =
      Test_n,

    Cases,

    Controls,

    Prevalence =
      round(
        Prevalence,
        3L
      ),

    AUROC =
      round(
        AUROC,
        3L
      ),

    `AUROC 95% CI` =
      sprintf(
        "%.3f\u2013%.3f",
        AUROC_lower95,
        AUROC_upper95
      ),

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
      )
  )
]


publication_loco_table[
  ,
  Cohort := factor(
    Cohort,
    levels = cohort_order
  )
]


data.table::setorder(
  publication_loco_table,
  Cohort
)


print(
  publication_loco_table
)


data.table::fwrite(
  publication_loco_table,
  file.path(
    TABLE_DIR,
    "Table_2_locked_model_strict_LOCO_performance.csv"
  )
)
