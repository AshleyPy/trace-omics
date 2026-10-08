# ======================================================================
# PART 7
# FIGURE 4
# STRICT LOCO CALIBRATION
#
# POINTS:
#   empirical observed case proportions
#
# LINE:
#   empirical binned calibration curve
#
# RIBBON:
#   exact binomial 95% CI for observed case proportion
# ======================================================================


calibration_rows <- list()


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


  calibration_rows[[cohort_index]] <-
    make_calibration_data(

      truth =
        current_dt$truth,

      probability =
        current_dt$probability_case,

      label =
        current_cohort,

      number_of_bins =
        6L
    )
}


calibration_dt <- data.table::rbindlist(
  calibration_rows,
  use.names = TRUE,
  fill = TRUE
)


calibration_dt[
  ,
  Group := factor(
    Group,
    levels = cohort_order
  )
]


data.table::setorder(
  calibration_dt,
  Group,
  Mean_predicted
)


# ======================================================================
# FIGURE
# ======================================================================

p_loco_calibration <- ggplot2::ggplot() +

  # --------------------------------------------------------------------
  # Exact binomial 95% CI ribbons
  # --------------------------------------------------------------------

  ggplot2::geom_ribbon(

    data =
      calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      ymin = Lower95,
      ymax = Upper95,
      fill = Group,
      group = Group
    ),

    alpha =
      0.14,

    colour =
      NA,

    inherit.aes =
      FALSE
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
  # Empirical calibration curve
  # --------------------------------------------------------------------

  ggplot2::geom_line(

    data =
      calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      y = Observed,
      colour = Group,
      group = Group
    ),

    linewidth =
      1.00
  ) +

  ggplot2::geom_point(

    data =
      calibration_dt,

    ggplot2::aes(
      x = Mean_predicted,
      y = Observed,
      colour = Group
    ),

    size =
      2.8
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

  ggplot2::coord_equal() +

  ggplot2::labs(

    title =
      "Calibration in completely unseen TB cohorts",

    subtitle =
      "Binned observed probabilities with exact binomial 95% confidence bands",

    x =
      "Mean predicted probability",

    y =
      "Observed case proportion"
  ) +

  publication_theme()


print(
  p_loco_calibration
)


# ======================================================================
# SAVE
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      FIGURE_DIR,
      "Figure_4_strict_LOCO_calibration_95CI.pdf"
    ),

  plot =
    p_loco_calibration,

  width =
    7.6,

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
      "Figure_4_strict_LOCO_calibration_95CI.png"
    ),

  plot =
    p_loco_calibration,

  width =
    7.6,

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
  calibration_dt,
  file.path(
    FIGURE_DIR,
    "Figure_4_calibration_plot_data.csv"
  )
)
