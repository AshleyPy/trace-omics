# ======================================================================
# PART 19A
# LATEX TABLE 1
# ======================================================================

writeLines(

  capture.output(

    knitr::kable(

      publication_optimisation_table,

      format = "latex",

      booktabs = TRUE,

      caption = paste0(
        "Development-only optimisation of the rescued TB classifier. ",
        "Candidate configurations were ranked primarily by worst-cohort ",
        "leave-one-cohort-out AUROC."
      ),

      label = "tab:tb_rescue_optimisation"
    )
  ),

  file.path(
    TABLE_DIR,
    "Table_1_model_optimisation.tex"
  )
)


# ======================================================================
# PART 19B
# LATEX TABLE 2
# ======================================================================

writeLines(

  capture.output(

    knitr::kable(

      publication_loco_table,

      format = "latex",

      booktabs = TRUE,

      caption = paste0(
        "Strict leave-one-cohort-out performance of the locked TB classifier. ",
        "Each target cohort was excluded completely from feature selection ",
        "and classifier fitting."
      ),

      label = "tab:tb_strict_loco"
    )
  ),

  file.path(
    TABLE_DIR,
    "Table_2_locked_model_strict_LOCO_performance.tex"
  )
)


# ======================================================================
# PART 19C
# LATEX TABLE 3
# ======================================================================

table3_latex <- diagnostic_ice_publication_result[
  ,
  .(

    Representation,

    Genes =
      TOP_N,

    Model,

    `Class weighting` =
      ifelse(
        Class_weighting,
        "Yes",
        "No"
      ),

    Threshold =
      round(
        Threshold,
        2L
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


writeLines(

  capture.output(

    knitr::kable(

      table3_latex,

      format = "latex",

      booktabs = TRUE,

      caption = paste0(
        "Diagnostic held-out evaluation of the rescued TB classifier. ",
        "This population is reported as diagnostic rather than pristine ",
        "confirmatory validation because model rescue followed inspection ",
        "of the original test failure."
      ),

      label = "tab:tb_diagnostic_ice"
    )
  ),

  file.path(
    TABLE_DIR,
    "Table_3_diagnostic_ICE_performance.tex"
  )
)
