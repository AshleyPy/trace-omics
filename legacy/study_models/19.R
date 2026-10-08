# ======================================================================
# PART 17
# ALE FEATURE-EFFECT PROFILES
# ======================================================================

profile_genes <- importance_summary_dt[
  seq_len(
    PROFILE_N
  ),
  Gene
]


cat(
  "\nGenes selected for ALE profiles:\n"
)

print(
  profile_genes
)


set.seed(
  SEED_PUBLICATION
)


dalex_ale <- tryCatch(

  ingredients::accumulated_dependency(
    tb_rf_explainer,
    variables = profile_genes,
    N = NULL
  ),

  error = function(e) {

    warning(
      "ALE calculation failed: ",
      conditionMessage(
        e
      )
    )

    NULL
  }
)


if (!is.null(dalex_ale)) {

  saveRDS(
    dalex_ale,
    file.path(
      DALEX_DIR,
      "TB_final_RF_top_gene_ALE.rds"
    )
  )

  data.table::fwrite(
    data.table::as.data.table(
      dalex_ale
    ),
    file.path(
      DALEX_DIR,
      "TB_final_RF_top_gene_ALE.csv"
    )
  )

  p_dalex_ale <- plot(
    dalex_ale
  ) +

    ggplot2::labs(
      title = "Expression-response profiles of the most predictive TB genes",
      subtitle = "Accumulated local effects from the locked Random Forest"
    ) +

    publication_theme()

  print(
    p_dalex_ale
  )

  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR,
      "Figure_7_top_gene_ALE_profiles.pdf"
    ),
    plot = p_dalex_ale,
    width = 9.0,
    height = 6.5,
    units = "in",
    dpi = 600,
    device = cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(
      FIGURE_DIR,
      "Figure_7_top_gene_ALE_profiles.png"
    ),
    plot = p_dalex_ale,
    width = 9.0,
    height = 6.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}
