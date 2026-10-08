# ======================================================================
# PART 16
# FIGURE 6
# PREDICTIVE GENE IMPORTANCE
# ======================================================================

importance_plot_dt <- data.table::copy(
  importance_summary_dt[
    1:min(
      TOP_IMPORTANCE_GENES,
      .N
    )
  ]
)


importance_plot_dt[
  ,
  Gene := factor(
    Gene,
    levels = rev(
      Gene
    )
  )
]


p_dalex_importance <- ggplot2::ggplot(
  importance_plot_dt,
  ggplot2::aes(
    x = Gene,
    y = Mean_Delta_AUROC
  )
) +

  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.55
  ) +

  ggplot2::geom_col(
    width = 0.70,
    fill = SET1[[1L]],
    colour = "black",
    linewidth = 0.25
  ) +

  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = Lower95_Delta_AUROC,
      ymax = Upper95_Delta_AUROC
    ),
    width = 0.16,
    linewidth = 0.40
  ) +

  ggplot2::coord_flip() +

  ggplot2::labs(
    title = "Predictive genes used by the locked TB classifier",
    subtitle = paste0(
      "Post-hoc permutation importance; ",
      B_IMPORTANCE,
      " permutations per gene"
    ),
    x = NULL,
    y = expression(
      "Permutation importance (" *
        Delta *
        "AUROC)"
    )
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(
      face = "italic",
      size = 9
    )
  )


print(
  p_dalex_importance
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_6_predictive_gene_importance.pdf"
  ),
  plot = p_dalex_importance,
  width = 7.3,
  height = 6.3,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_6_predictive_gene_importance.png"
  ),
  plot = p_dalex_importance,
  width = 7.3,
  height = 6.3,
  units = "in",
  dpi = 600,
  bg = "white"
)
