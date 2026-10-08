# ======================================================================
# PART 10A
# DIAGNOSTIC ICE PR
# EMPIRICAL PR + STRATIFIED-BOOTSTRAP 95% CI
# ======================================================================


ice_pr_object <- PRROC::pr.curve(

  scores.class0 =
    diagnostic_ice_predictions$probability_case[
      diagnostic_ice_predictions$truth ==
        POSITIVE_CLASS
    ],

  scores.class1 =
    diagnostic_ice_predictions$probability_case[
      diagnostic_ice_predictions$truth ==
        NEGATIVE_CLASS
    ],

  curve =
    TRUE
)


ice_prevalence <- mean(
  diagnostic_ice_predictions$truth ==
    POSITIVE_CLASS
)


ice_pr_dt <- data.table::data.table(

  Recall =
    as.numeric(
      ice_pr_object$curve[
        ,
        1L
      ]
    ),

  Precision =
    as.numeric(
      ice_pr_object$curve[
        ,
        2L
      ]
    )
)


ice_pr_dt <- ice_pr_dt[
  is.finite(
    Recall
  ) &
    is.finite(
      Precision
    )
]


# ======================================================================
# ICE PR 95% CI
# ======================================================================

ice_pr_ribbon_dt <- bootstrap_pr_ci(

  truth =
    diagnostic_ice_predictions$truth,

  probability =
    diagnostic_ice_predictions$probability_case,

  recall_grid =
    PR_RECALL_GRID,

  bootstraps =
    PR_BOOTSTRAPS,

  seed =
    SEED_PUBLICATION +
      90000L
)


ice_pr_ribbon_dt[
  ,
  `:=`(

    Precision_lower =
      pmax(
        0,
        pmin(
          1,
          Precision_lower
        )
      ),

    Precision_upper =
      pmax(
        0,
        pmin(
          1,
          Precision_upper
        )
      )
  )
]


p_ice_pr <- ggplot2::ggplot() +

  # --------------------------------------------------------------------
  # 95% CI
  # --------------------------------------------------------------------

  ggplot2::geom_ribbon(

    data =
      ice_pr_ribbon_dt,

    ggplot2::aes(
      x = Recall,
      ymin = Precision_lower,
      ymax = Precision_upper
    ),

    fill =
      SET1[[1L]],

    alpha =
      0.15,

    colour =
      NA
  ) +

  # --------------------------------------------------------------------
  # Prevalence
  # --------------------------------------------------------------------

  ggplot2::geom_hline(

    yintercept =
      ice_prevalence,

    linetype =
      "dashed",

    colour =
      SET1[[1L]],

    linewidth =
      0.75
  ) +

  # --------------------------------------------------------------------
  # Empirical PR
  # --------------------------------------------------------------------

  ggplot2::geom_path(

    data =
      ice_pr_dt,

    ggplot2::aes(
      x = Recall,
      y = Precision
    ),

    linewidth =
      1.10,

    colour =
      SET1[[1L]]
  ) +

  ggplot2::scale_x_continuous(

    limits = c(
      0,
      1
    ),

    breaks = seq(
      0,
      1,
      by = 0.20
    ),

    expand = c(
      0,
      0
    )
  ) +

  ggplot2::scale_y_continuous(

    limits = c(
      0,
      1
    ),

    breaks = seq(
      0,
      1,
      by = 0.20
    ),

    expand = c(
      0,
      0
    )
  ) +

  ggplot2::labs(

    title =
      "Diagnostic held-out precision-recall",

    subtitle =
      sprintf(
        "AUPRC %.3f | prevalence %.3f | pointwise bootstrap 95%% CI",
        ice_pr_object$auc.integral,
        ice_prevalence
      ),

    x =
      "Recall",

    y =
      "Precision"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none"
  )


# ======================================================================
# PART 10B
# DIAGNOSTIC ICE CALIBRATION + 95% CI RIBBON
# ======================================================================


ice_calibration_dt <- make_calibration_data(

  truth =
    diagnostic_ice_predictions$truth,

  probability =
    diagnostic_ice_predictions$probability_case,

  label =
    "Diagnostic ICE",

  number_of_bins =
    6L
)


data.table::setorder(
  ice_calibration_dt,
  Mean_predicted
)


p_ice_calibration <- ggplot2::ggplot() +

  # --------------------------------------------------------------------
  # Exact binomial 95% CI
  # --------------------------------------------------------------------

  ggplot2::geom_ribbon(

    data =
      ice_calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      ymin = Lower95,
      ymax = Upper95
    ),

    fill =
      SET1[[1L]],

    alpha =
      0.15,

    colour =
      NA
  ) +

  # --------------------------------------------------------------------
  # Perfect calibration
  # --------------------------------------------------------------------

  ggplot2::geom_abline(

    intercept =
      0,

    slope =
      1,

    linetype =
      "dashed",

    colour =
      "grey45",

    linewidth =
      0.70
  ) +

  # --------------------------------------------------------------------
  # Empirical calibration
  # --------------------------------------------------------------------

  ggplot2::geom_line(

    data =
      ice_calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      y = Observed
    ),

    linewidth =
      0.95,

    colour =
      SET1[[1L]]
  ) +

  ggplot2::geom_point(

    data =
      ice_calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      y = Observed
    ),

    size =
      2.8,

    colour =
      SET1[[1L]]
  ) +

  ggplot2::scale_x_continuous(

    limits = c(
      0,
      1
    ),

    breaks = seq(
      0,
      1,
      by = 0.20
    ),

    expand = c(
      0,
      0
    )
  ) +

  ggplot2::scale_y_continuous(

    limits = c(
      0,
      1
    ),

    breaks = seq(
      0,
      1,
      by = 0.20
    ),

    expand = c(
      0,
      0
    )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(

    title =
      "Diagnostic held-out calibration",

    subtitle =
      "Binned observed case proportions with exact binomial 95% confidence band",

    x =
      "Mean predicted probability",

    y =
      "Observed case proportion"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none"
  )


# ======================================================================
# PART 10C
# COMBINED DIAGNOSTIC ICE PANEL
#
# p_ice_roc comes from Part 9.
# ======================================================================


p_ice_panel <- (

  p_ice_roc +

    p_ice_pr +

    p_ice_calibration

) +

  patchwork::plot_layout(
    ncol = 3
  ) +

  patchwork::plot_annotation(

    title =
      "Diagnostic evaluation of the rescued TB classifier",

    tag_levels =
      "A"
  )


print(
  p_ice_panel
)


# ======================================================================
# SAVE
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      FIGURE_DIR,
      "Figure_S1_diagnostic_ICE_ROC_PR_calibration_95CI.pdf"
    ),

  plot =
    p_ice_panel,

  width =
    14.8,

  height =
    5.2,

  units =
    "in",

  dpi =
    600,

  device =
    cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      FIGURE_DIR,
      "Figure_S1_diagnostic_ICE_ROC_PR_calibration_95CI.png"
    ),

  plot =
    p_ice_panel,

  width =
    14.8,

  height =
    5.2,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


data.table::fwrite(
  ice_pr_ribbon_dt,
  file.path(
    FIGURE_DIR,
    "Figure_S1_ICE_PR_95CI_ribbon_data.csv"
  )
)


data.table::fwrite(
  ice_calibration_dt,
  file.path(
    FIGURE_DIR,
    "Figure_S1_ICE_calibration_plot_data.csv"
  )
)
