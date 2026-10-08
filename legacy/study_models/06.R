# ======================================================================
# PART 4
# FIGURE 1
# DEVELOPMENT-ONLY OPTIMISATION
# ======================================================================

figure1_dt <- data.table::copy(
  optimisation_metrics
)

figure1_dt[
  ,
  Feature_set := factor(
    paste0(
      TOP_N,
      " genes"
    ),
    levels = c(
      "25 genes",
      "50 genes",
      "100 genes"
    )
  )
]

figure1_dt[
  ,
  Weighting := factor(
    ifelse(
      Class_weighting,
      "Class weighted",
      "Unweighted"
    ),
    levels = c(
      "Unweighted",
      "Class weighted"
    )
  )
]

figure1_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = cohort_order
  )
]


p_optimisation <- ggplot2::ggplot(
  figure1_dt,
  ggplot2::aes(
    x = Target_cohort,
    y = AUROC,
    group = interaction(
      Feature_set,
      Weighting
    ),
    colour = Feature_set,
    shape = Weighting
  )
) +

  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.60
  ) +

  ggplot2::geom_line(
    linewidth = 0.70,
    alpha = 0.80
  ) +

  ggplot2::geom_point(
    size = 3.2
  ) +

  ggplot2::scale_colour_manual(
    values = SET1[
      1:3
    ]
  ) +

  ggplot2::coord_cartesian(
    ylim = c(
      0.45,
      1.00
    )
  ) +

  ggplot2::labs(
    title = "Development-only optimisation of the TB classifier",
    subtitle = "Model selection prioritised worst-cohort strict leave-one-cohort-out AUROC",
    x = "Completely held-out cohort",
    y = "AUROC",
    colour = "Feature set",
    shape = "Classifier weighting"
  ) +

  publication_theme()


print(
  p_optimisation
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_1_development_optimisation.pdf"
  ),
  plot = p_optimisation,
  width = 7.8,
  height = 5.3,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_1_development_optimisation.png"
  ),
  plot = p_optimisation,
  width = 7.8,
  height = 5.3,
  units = "in",
  dpi = 600,
  bg = "white"
)
