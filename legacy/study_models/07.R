# ======================================================================
# PART 5
# FIGURE 2
# STRICT LOCO EMPIRICAL STEPWISE ROC + 95% CI RIBBONS
#
# ROC LINE:
#   empirical ROC coordinates
#   geom_step(direction = "vh")
#
# AUROC:
#   empirical
#
# AUROC 95% CI:
#   DeLong
#
# RIBBON:
#   pointwise bootstrap 95% sensitivity interval
#   at fixed specificity
#
# IMPORTANT:
#   NO ROC SMOOTHING IS USED.
#
# No paired DeLong test is performed between ACS, Catalysis and GC6
# because these curves arise from different observations.
# ======================================================================


roc_curve_rows <- list()

roc_ribbon_rows <- list()

roc_summary_rows <- list()

roc_legend_labels <- character()


for (
  cohort_index in seq_along(
    cohort_order
  )
) {

  current_cohort <- cohort_order[[cohort_index]]

  current_dt <- locked_loco_predictions[
    Target_cohort == current_cohort
  ]

  current_truth <- factor(
    current_dt$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )


  # --------------------------------------------------------------------
  # Empirical ROC object
  # --------------------------------------------------------------------

  current_roc <- pROC::roc(
    response = current_truth,
    predictor = current_dt$probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )


  # --------------------------------------------------------------------
  # Empirical AUROC
  # --------------------------------------------------------------------

  current_auc <- as.numeric(
    pROC::auc(
      current_roc
    )
  )


  # --------------------------------------------------------------------
  # DeLong 95% AUROC CI
  # --------------------------------------------------------------------

  current_auc_ci <- as.numeric(
    pROC::ci.auc(
      current_roc,
      method = "delong",
      conf.level = 0.95
    )
  )


  # --------------------------------------------------------------------
  # Empirical ROC coordinates
  # --------------------------------------------------------------------

  current_coordinates <- pROC::coords(
    current_roc,
    x = "all",
    ret = c(
      "specificity",
      "sensitivity"
    ),
    transpose = FALSE
  )

  current_curve_dt <- data.table::data.table(

    Target_cohort = current_cohort,

    False_positive_rate =
      1 -
      as.numeric(
        current_coordinates$specificity
      ),

    True_positive_rate =
      as.numeric(
        current_coordinates$sensitivity
      )
  )

  current_curve_dt <- current_curve_dt[
    is.finite(
      False_positive_rate
    ) &
      is.finite(
        True_positive_rate
      )
  ]

  data.table::setorder(
    current_curve_dt,
    False_positive_rate,
    True_positive_rate
  )

  roc_curve_rows[[cohort_index]] <-
    current_curve_dt


  # --------------------------------------------------------------------
  # Pointwise bootstrap 95% ROC CI
  # --------------------------------------------------------------------

  set.seed(
    SEED_PUBLICATION +
      cohort_index *
        1000L
  )

  current_curve_ci <- pROC::ci.se(
    current_roc,
    specificities = ROC_SPECIFICITY_GRID,
    conf.level = 0.95,
    boot.n = ROC_BOOTSTRAPS,
    boot.stratified = TRUE
  )

  current_ribbon_dt <- data.table::data.table(

    Target_cohort = current_cohort,

    False_positive_rate =
      1 -
      ROC_SPECIFICITY_GRID,

    Sensitivity_lower =
      as.numeric(
        current_curve_ci[
          ,
          1L
        ]
      ),

    Sensitivity_median =
      as.numeric(
        current_curve_ci[
          ,
          2L
        ]
      ),

    Sensitivity_upper =
      as.numeric(
        current_curve_ci[
          ,
          3L
        ]
      )
  )

  current_ribbon_dt[
    ,
    Sensitivity_lower := pmax(
      0,
      pmin(
        1,
        Sensitivity_lower
      )
    )
  ]

  current_ribbon_dt[
    ,
    Sensitivity_upper := pmax(
      0,
      pmin(
        1,
        Sensitivity_upper
      )
    )
  ]

  data.table::setorder(
    current_ribbon_dt,
    False_positive_rate
  )

  roc_ribbon_rows[[cohort_index]] <-
    current_ribbon_dt


  # --------------------------------------------------------------------
  # Legend
  # --------------------------------------------------------------------

  current_label <- paste0(
    current_cohort,
    ": ",
    sprintf(
      "%.3f",
      current_auc
    ),
    " (",
    sprintf(
      "%.3f",
      current_auc_ci[[1L]]
    ),
    "\u2013",
    sprintf(
      "%.3f",
      current_auc_ci[[3L]]
    ),
    ")"
  )

  roc_legend_labels[[current_cohort]] <-
    current_label


  # --------------------------------------------------------------------
  # Summary
  # --------------------------------------------------------------------

  roc_summary_rows[[cohort_index]] <-
    data.table::data.table(

      Cohort = current_cohort,

      n = nrow(
        current_dt
      ),

      Controls = sum(
        current_truth ==
          NEGATIVE_CLASS
      ),

      Cases = sum(
        current_truth ==
          POSITIVE_CLASS
      ),

      AUROC = current_auc,

      AUROC_lower95 =
        current_auc_ci[[1L]],

      AUROC_upper95 =
        current_auc_ci[[3L]]
    )
}


# ======================================================================
# COMBINE
# ======================================================================

roc_curve_dt <- data.table::rbindlist(
  roc_curve_rows,
  use.names = TRUE,
  fill = TRUE
)

roc_ribbon_dt <- data.table::rbindlist(
  roc_ribbon_rows,
  use.names = TRUE,
  fill = TRUE
)

roc_summary <- data.table::rbindlist(
  roc_summary_rows,
  use.names = TRUE,
  fill = TRUE
)


roc_curve_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = cohort_order
  )
]

roc_ribbon_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = cohort_order
  )
]


data.table::setorder(
  roc_curve_dt,
  Target_cohort,
  False_positive_rate,
  True_positive_rate
)

data.table::setorder(
  roc_ribbon_dt,
  Target_cohort,
  False_positive_rate
)


# ======================================================================
# PLOT
# ======================================================================

p_loco_roc <- ggplot2::ggplot() +

  # --------------------------------------------------------------------
  # Pointwise 95% bootstrap confidence region
  # --------------------------------------------------------------------

  ggplot2::geom_ribbon(
    data = roc_ribbon_dt,
    ggplot2::aes(
      x = False_positive_rate,
      ymin = Sensitivity_lower,
      ymax = Sensitivity_upper,
      fill = Target_cohort,
      group = Target_cohort
    ),
    alpha = 0.14,
    colour = NA,
    inherit.aes = FALSE
  ) +

  # --------------------------------------------------------------------
  # Chance
  # --------------------------------------------------------------------

  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.70
  ) +

  # --------------------------------------------------------------------
  # EMPIRICAL STEPWISE ROC
  # --------------------------------------------------------------------

  ggplot2::geom_step(
    data = roc_curve_dt,
    ggplot2::aes(
      x = False_positive_rate,
      y = True_positive_rate,
      colour = Target_cohort,
      group = Target_cohort
    ),
    linewidth = 1.15,
    direction = "vh"
  ) +

  # --------------------------------------------------------------------
  # Colours
  # --------------------------------------------------------------------

  ggplot2::scale_colour_manual(
    values = cohort_colours,
    limits = cohort_order,
    breaks = cohort_order,
    labels = roc_legend_labels[
      cohort_order
    ],
    drop = FALSE,
    name = "Held-out cohort \u2014 AUROC (95% CI)"
  ) +

  ggplot2::scale_fill_manual(
    values = cohort_colours,
    limits = cohort_order,
    breaks = cohort_order,
    drop = FALSE,
    guide = "none"
  ) +

  # --------------------------------------------------------------------
  # Axes
  # --------------------------------------------------------------------

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
    labels = scales::label_number(
      accuracy = 0.1
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
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = c(
      0,
      0
    )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(
    title = "Strict cross-cohort discrimination of the locked TB classifier",
    subtitle = "Empirical ROC curves with pointwise 95% bootstrap confidence bands",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +

  publication_theme()


print(
  p_loco_roc
)


# ======================================================================
# SAVE
# ======================================================================

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_2_strict_LOCO_ROC_STEPWISE_95CI.pdf"
  ),
  plot = p_loco_roc,
  width = 8.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_2_strict_LOCO_ROC_STEPWISE_95CI.png"
  ),
  plot = p_loco_roc,
  width = 8.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

data.table::fwrite(
  roc_summary,
  file.path(
    TABLE_DIR,
    "Figure_2_ROC_summary.csv"
  )
)

data.table::fwrite(
  roc_curve_dt,
  file.path(
    FIGURE_DIR,
    "Figure_2_ROC_empirical_coordinates.csv"
  )
)

data.table::fwrite(
  roc_ribbon_dt,
  file.path(
    FIGURE_DIR,
    "Figure_2_ROC_95CI_ribbon_data.csv"
  )
)
