# ======================================================================
# PART 6
# FIGURE 3
# STRICT LOCO PRECISION-RECALL
#
# CURVE:
#   empirical PR curve
#
# RIBBON:
#   pointwise stratified-bootstrap 95% CI
#
# HORIZONTAL DASHED LINE:
#   held-out cohort case prevalence
#   shown using the same cohort colour
# ======================================================================


pr_curve_rows <- list()

pr_ribbon_rows <- list()

pr_summary_rows <- list()

pr_legend_labels <- character()


for (
  cohort_index in seq_along(
    cohort_order
  )
) {

  current_cohort <- cohort_order[[cohort_index]]


  current_dt <- locked_loco_predictions[
    Target_cohort ==
      current_cohort
  ]


  current_truth <- as.character(
    current_dt$truth
  )


  # ====================================================================
  # EMPIRICAL PR CURVE
  # ====================================================================

  current_pr <- PRROC::pr.curve(

    scores.class0 =
      current_dt$probability_case[
        current_truth ==
          POSITIVE_CLASS
      ],

    scores.class1 =
      current_dt$probability_case[
        current_truth ==
          NEGATIVE_CLASS
      ],

    curve =
      TRUE
  )


  current_prevalence <- mean(
    current_truth ==
      POSITIVE_CLASS
  )


  current_curve_dt <- data.table::data.table(

    Target_cohort =
      current_cohort,

    Recall =
      as.numeric(
        current_pr$curve[
          ,
          1L
        ]
      ),

    Precision =
      as.numeric(
        current_pr$curve[
          ,
          2L
        ]
      )
  )


  current_curve_dt <- current_curve_dt[
    is.finite(
      Recall
    ) &
      is.finite(
        Precision
      )
  ]


  pr_curve_rows[[cohort_index]] <-
    current_curve_dt


  # ====================================================================
  # STRATIFIED BOOTSTRAP 95% CI RIBBON
  # ====================================================================

  current_ribbon_dt <- bootstrap_pr_ci(

    truth =
      current_truth,

    probability =
      current_dt$probability_case,

    recall_grid =
      PR_RECALL_GRID,

    bootstraps =
      PR_BOOTSTRAPS,

    seed =
      SEED_PUBLICATION +
        cohort_index *
          10000L
  )


  current_ribbon_dt[
    ,
    Target_cohort :=
      current_cohort
  ]


  current_ribbon_dt[
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


  pr_ribbon_rows[[cohort_index]] <-
    current_ribbon_dt


  # ====================================================================
  # LEGEND LABEL
  # ====================================================================

  current_label <- paste0(

    current_cohort,

    ": AUPRC ",

    sprintf(
      "%.3f",
      current_pr$auc.integral
    ),

    " | prevalence ",

    sprintf(
      "%.3f",
      current_prevalence
    )
  )


  pr_legend_labels[[current_cohort]] <-
    current_label


  # ====================================================================
  # SUMMARY
  # ====================================================================

  pr_summary_rows[[cohort_index]] <-
    data.table::data.table(

      Cohort =
        current_cohort,

      n =
        nrow(
          current_dt
        ),

      Cases =
        sum(
          current_truth ==
            POSITIVE_CLASS
        ),

      Controls =
        sum(
          current_truth ==
            NEGATIVE_CLASS
        ),

      Prevalence =
        current_prevalence,

      AUPRC =
        current_pr$auc.integral
    )
}


# ======================================================================
# COMBINE
# ======================================================================

pr_curve_dt <- data.table::rbindlist(
  pr_curve_rows,
  use.names = TRUE,
  fill = TRUE
)


pr_ribbon_dt <- data.table::rbindlist(
  pr_ribbon_rows,
  use.names = TRUE,
  fill = TRUE
)


pr_summary <- data.table::rbindlist(
  pr_summary_rows,
  use.names = TRUE,
  fill = TRUE
)


pr_curve_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = cohort_order
  )
]


pr_ribbon_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = cohort_order
  )
]


pr_summary[
  ,
  Cohort := factor(
    Cohort,
    levels = cohort_order
  )
]


data.table::setorder(
  pr_ribbon_dt,
  Target_cohort,
  Recall
)


# ======================================================================
# FIGURE
# ======================================================================

p_loco_pr <- ggplot2::ggplot() +

  # --------------------------------------------------------------------
  # 95% PR confidence ribbons
  # --------------------------------------------------------------------

  ggplot2::geom_ribbon(

    data =
      pr_ribbon_dt,

    ggplot2::aes(
      x = Recall,
      ymin = Precision_lower,
      ymax = Precision_upper,
      fill = Target_cohort,
      group = Target_cohort
    ),

    alpha =
      0.14,

    colour =
      NA,

    inherit.aes =
      FALSE
  ) +

  # --------------------------------------------------------------------
  # COHORT-SPECIFIC PREVALENCE LINES
  # --------------------------------------------------------------------

  ggplot2::geom_hline(

    data =
      pr_summary,

    ggplot2::aes(
      yintercept = Prevalence,
      colour = Cohort
    ),

    linetype =
      "dashed",

    linewidth =
      0.70,

    alpha =
      0.90,

    show.legend =
      FALSE
  ) +

  # --------------------------------------------------------------------
  # Empirical PR curves
  # --------------------------------------------------------------------

  ggplot2::geom_path(

    data =
      pr_curve_dt,

    ggplot2::aes(
      x = Recall,
      y = Precision,
      colour = Target_cohort,
      group = Target_cohort
    ),

    linewidth =
      1.10
  ) +

  # --------------------------------------------------------------------
  # Colours
  # --------------------------------------------------------------------

  ggplot2::scale_colour_manual(

    values =
      cohort_colours,

    limits =
      cohort_order,

    breaks =
      cohort_order,

    labels =
      pr_legend_labels[
        cohort_order
      ],

    drop =
      FALSE,

    name =
      "Held-out cohort"
  ) +

  ggplot2::scale_fill_manual(

    values =
      cohort_colours,

    limits =
      cohort_order,

    breaks =
      cohort_order,

    guide =
      "none"
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

    labels =
      scales::label_number(
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

    labels =
      scales::label_number(
        accuracy = 0.1
      ),

    expand = c(
      0,
      0
    )
  ) +

  ggplot2::labs(

    title =
      "Precision-recall performance in completely unseen TB cohorts",

    subtitle =
      paste0(
        "Empirical PR curves with pointwise 95% stratified-bootstrap confidence bands; ",
        "dashed lines show cohort-specific prevalence"
      ),

    x =
      "Recall",

    y =
      "Precision"
  ) +

  publication_theme()


print(
  p_loco_pr
)


# ======================================================================
# SAVE
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      FIGURE_DIR,
      "Figure_3_strict_LOCO_PR_95CI.pdf"
    ),

  plot =
    p_loco_pr,

  width =
    8.3,

  height =
    5.8,

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
      "Figure_3_strict_LOCO_PR_95CI.png"
    ),

  plot =
    p_loco_pr,

  width =
    8.3,

  height =
    5.8,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


data.table::fwrite(
  pr_summary,
  file.path(
    TABLE_DIR,
    "Figure_3_PR_summary.csv"
  )
)


data.table::fwrite(
  pr_ribbon_dt,
  file.path(
    FIGURE_DIR,
    "Figure_3_PR_95CI_ribbon_data.csv"
  )
)
