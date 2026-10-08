# ======================================================================
# PART 8
# FIGURE 5
# MAIN STRICT-LOCO VALIDATION PANEL
#
# A = ROC + 95% CI
# B = PR + 95% CI + prevalence
# C = Calibration + 95% CI
# ======================================================================


p_main_validation <- (

  p_loco_roc +

    p_loco_pr +

    p_loco_calibration

) +

  patchwork::plot_layout(

    ncol =
      3,

    widths = c(
      1.08,
      1.08,
      1.00
    ),

    guides =
      "keep"
  ) +

  patchwork::plot_annotation(

    title =
      "Transportability and predictive reliability of the locked TB transcriptomic classifier",

    tag_levels =
      "A"
  )


print(
  p_main_validation
)


ggplot2::ggsave(

  filename =
    file.path(
      FIGURE_DIR,
      "Figure_5_MAIN_LOCO_ROC_PR_calibration_95CI.pdf"
    ),

  plot =
    p_main_validation,

  width =
    16.5,

  height =
    5.6,

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
      "Figure_5_MAIN_LOCO_ROC_PR_calibration_95CI.png"
    ),

  plot =
    p_main_validation,

  width =
    16.5,

  height =
    5.6,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)
