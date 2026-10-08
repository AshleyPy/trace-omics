# ======================================================================
# RA LIGHTGBM vs TABPFN
# FINAL DALEX IMPORTANCE + ALE FIGURES
#
# START HERE AFTER BOTH PERMUTATION IMPORTANCE RUNS ARE COMPLETE
#
# DOES NOT RERUN:
#   - diagnostic ICE
#   - model fitting
#   - feature selection
#   - tuning
#   - LightGBM permutation importance
#   - TabPFN permutation importance
# ======================================================================


# ======================================================================
# 1. PACKAGES
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(ingredients)
  library(DALEX)
})


# ======================================================================
# 2. PATHS
# ======================================================================

PROJECT_DIR <- file.path(
  Sys.getenv("USERPROFILE"),
  "Desktop",
  "RA_FULL_DATA_FINAL_MODEL"
)


FOUNDATION_DIR <- file.path(
  PROJECT_DIR,
  "103_LIGHTGBM_VS_TABPFN"
)


DALEX_COMPARE_DIR <- file.path(
  FOUNDATION_DIR,
  "DALEX_LightGBM_vs_TabPFN"
)


# ======================================================================
# 3. LOAD COMPLETED IMPORTANCE RESULTS
# ======================================================================

lgbm_importance_file <- file.path(
  DALEX_COMPARE_DIR,
  "RA_LightGBM_permutation_importance_summary.csv"
)


tabpfn_importance_file <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_summary.csv"
)


if (!file.exists(lgbm_importance_file)) {
  stop(
    "LightGBM importance summary not found:\n",
    lgbm_importance_file
  )
}


if (!file.exists(tabpfn_importance_file)) {
  stop(
    "TabPFN importance summary not found:\n",
    tabpfn_importance_file
  )
}


lgbm_importance_summary <- data.table::fread(
  lgbm_importance_file
)


tabpfn_importance_summary <- data.table::fread(
  tabpfn_importance_file
)


lgbm_importance_summary[
  ,
  Model := "LightGBM"
]


tabpfn_importance_summary[
  ,
  Model := "TabPFN"
]


data.table::setorder(
  lgbm_importance_summary,
  Rank
)


data.table::setorder(
  tabpfn_importance_summary,
  Rank
)


# ======================================================================
# 4. HARD AUDIT
# ======================================================================

stopifnot(
  nrow(lgbm_importance_summary) == 100L
)


stopifnot(
  nrow(tabpfn_importance_summary) == 100L
)


stopifnot(
  data.table::uniqueN(lgbm_importance_summary$Gene) == 100L
)


stopifnot(
  data.table::uniqueN(tabpfn_importance_summary$Gene) == 100L
)


cat(
  "\n============================================================\n",
  "COMPLETED RA IMPORTANCE RESULTS LOADED\n",
  "============================================================\n",
  "LightGBM genes: ",
  nrow(lgbm_importance_summary),
  "\nTabPFN genes: ",
  nrow(tabpfn_importance_summary),
  "\nLightGBM baseline AUROC: ",
  sprintf(
    "%.3f",
    lgbm_importance_summary$Baseline_AUROC[[1L]]
  ),
  "\nTabPFN baseline AUROC: ",
  sprintf(
    "%.3f",
    tabpfn_importance_summary$Baseline_AUROC[[1L]]
  ),
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 5. COMBINE RESULTS
# ======================================================================

comparison_importance_summary <- data.table::rbindlist(
  list(
    lgbm_importance_summary,
    tabpfn_importance_summary
  ),
  use.names = TRUE,
  fill = TRUE
)


data.table::fwrite(
  comparison_importance_summary,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_vs_TabPFN_importance_combined.csv"
  )
)


# ======================================================================
# 6. AESTHETICS
# ======================================================================

SET1 <- RColorBrewer::brewer.pal(
  9,
  "Set1"
)


MODEL_ORDER <- c(
  "LightGBM",
  "TabPFN"
)


MODEL_COLOURS <- stats::setNames(
  SET1[1:2],
  MODEL_ORDER
)


publication_theme <- function() {

  ggplot2::theme_bw(
    base_size = 11
  ) +

    ggplot2::theme(

      plot.title =
        ggplot2::element_text(
          face = "bold",
          size = 13,
          colour = "black"
        ),

      plot.subtitle =
        ggplot2::element_text(
          size = 9.5,
          colour = "black"
        ),

      axis.title =
        ggplot2::element_text(
          face = "bold",
          colour = "black"
        ),

      axis.text =
        ggplot2::element_text(
          colour = "black"
        ),

      legend.title =
        ggplot2::element_text(
          face = "bold",
          colour = "black"
        ),

      legend.text =
        ggplot2::element_text(
          colour = "black"
        ),

      legend.position =
        "top",

      panel.grid.minor =
        ggplot2::element_blank(),

      strip.background =
        ggplot2::element_rect(
          fill = "white",
          colour = "black"
        ),

      strip.text =
        ggplot2::element_text(
          face = "bold",
          colour = "black"
        )
    )
}


# ======================================================================
# 7. TOP-20 IMPORTANCE PLOT FUNCTION
# ======================================================================

TOP_N_PLOT <- 20L


make_importance_plot <- function(
  importance_dt,
  model_name
) {

  plot_dt <- data.table::copy(
    importance_dt[
      Model == model_name
    ][
      1:min(
        TOP_N_PLOT,
        .N
      )
    ]
  )


  plot_dt[
    ,
    Gene := factor(
      Gene,
      levels = rev(Gene)
    )
  ]


  ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(
      x = Gene,
      y = Mean_Delta_AUROC,
      fill = Model
    )
  ) +

    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.5
    ) +

    ggplot2::geom_col(
      width = 0.70,
      colour = "black",
      linewidth = 0.2
    ) +

    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Lower95_Delta_AUROC,
        ymax = Upper95_Delta_AUROC
      ),
      width = 0.15,
      linewidth = 0.45
    ) +

    ggplot2::coord_flip() +

    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      name = NULL
    ) +

    ggplot2::labs(
      title = paste0(
        model_name,
        " feature importance"
      ),
      subtitle =
        "Repeated permutation importance in the diagnostic held-out cohort",
      x = NULL,
      y = expression(
        "Permutation importance (" *
          Delta *
          "AUROC)"
      )
    ) +

    publication_theme() +

    ggplot2::theme(
      axis.text.y =
        ggplot2::element_text(
          face = "italic",
          size = 9,
          colour = "black"
        ),
      panel.grid.major.y =
        ggplot2::element_blank()
    )
}


# ======================================================================
# 8. INDIVIDUAL IMPORTANCE FIGURES
# ======================================================================

p_lgbm_importance <- make_importance_plot(
  comparison_importance_summary,
  "LightGBM"
)


p_tabpfn_importance <- make_importance_plot(
  comparison_importance_summary,
  "TabPFN"
)


print(
  p_lgbm_importance
)


print(
  p_tabpfn_importance
)


# ======================================================================
# 9. COMBINED IMPORTANCE FIGURE
# ======================================================================

p_importance_combined <- (

  p_lgbm_importance +

    p_tabpfn_importance

) +

  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) +

  patchwork::plot_annotation(

    title =
      "Predictive feature reliance in rheumatoid arthritis",

    subtitle =
      paste0(
        "LightGBM and TabPFN evaluated in the identical ",
        "diagnostic held-out population"
      ),

    tag_levels =
      "A"
  ) &

  ggplot2::theme(
    legend.position = "top"
  )


print(
  p_importance_combined
)


# ======================================================================
# 10. SAVE IMPORTANCE FIGURES
# ======================================================================

ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_LightGBM_vs_TabPFN.pdf"
  ),
  plot = p_importance_combined,
  width = 13.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_LightGBM_vs_TabPFN.png"
  ),
  plot = p_importance_combined,
  width = 13.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 11. TOP-20 COMPARISON TABLE
# ======================================================================

top20_table <- data.table::rbindlist(
  list(

    lgbm_importance_summary[
      1:min(
        20L,
        .N
      )
    ],

    tabpfn_importance_summary[
      1:min(
        20L,
        .N
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)


data.table::fwrite(
  top20_table,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_top20_importance_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 12. CROSS-MODEL IMPORTANCE CONCORDANCE
# ======================================================================

importance_wide <- data.table::dcast(
  comparison_importance_summary,
  Gene ~ Model,
  value.var = "Mean_Delta_AUROC"
)


importance_spearman <- stats::cor(
  importance_wide$LightGBM,
  importance_wide$TabPFN,
  method = "spearman",
  use = "complete.obs"
)


lgbm_top20 <- lgbm_importance_summary[
  1:min(
    20L,
    .N
  ),
  Gene
]


tabpfn_top20 <- tabpfn_importance_summary[
  1:min(
    20L,
    .N
  ),
  Gene
]


shared_top20 <- intersect(
  lgbm_top20,
  tabpfn_top20
)


top20_union <- union(
  lgbm_top20,
  tabpfn_top20
)


top20_jaccard <- length(
  shared_top20
) /
  length(
    top20_union
  )


concordance_summary <- data.table::data.table(

  Genes_compared =
    nrow(
      importance_wide
    ),

  Spearman_rho =
    importance_spearman,

  Shared_top20 =
    length(
      shared_top20
    ),

  Top20_Jaccard =
    top20_jaccard
)


print(
  concordance_summary
)


data.table::fwrite(
  concordance_summary,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_importance_concordance.csv"
  )
)


data.table::fwrite(
  data.table::data.table(
    Gene = shared_top20
  ),
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_shared_top20_genes.csv"
  )
)


# ======================================================================
# 13. CONCORDANCE SCATTER
# ======================================================================

p_concordance <- ggplot2::ggplot(

  importance_wide,

  ggplot2::aes(
    x = LightGBM,
    y = TabPFN
  )
) +

  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey65"
  ) +

  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey65"
  ) +

  ggplot2::geom_point(
    size = 2.2,
    alpha = 0.75
  ) +

  ggplot2::geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.75
  ) +

  ggplot2::annotate(
    geom = "label",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.5,
    label = paste0(
      "Spearman \u03c1 = ",
      sprintf(
        "%.3f",
        importance_spearman
      ),
      "\nTop-20 Jaccard = ",
      sprintf(
        "%.3f",
        top20_jaccard
      )
    ),
    size = 3.4
  ) +

  ggplot2::labs(
    title =
      "Cross-model concordance of RA feature reliance",
    subtitle =
      "Each point represents one gene evaluated by both models",
    x = expression(
      "LightGBM importance (" *
        Delta *
        "AUROC)"
    ),
    y = expression(
      "TabPFN importance (" *
        Delta *
        "AUROC)"
    )
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none"
  )


print(
  p_concordance
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_importance_concordance.pdf"
  ),
  plot = p_concordance,
  width = 6.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_importance_concordance.png"
  ),
  plot = p_concordance,
  width = 6.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 14. SELECT TOP SIX GENES FOR ALE
# ======================================================================

PROFILE_N <- 6L


LGBM_PROFILE_GENES <- lgbm_importance_summary[
  1:min(
    PROFILE_N,
    .N
  ),
  Gene
]


TABPFN_PROFILE_GENES <- tabpfn_importance_summary[
  1:min(
    PROFILE_N,
    .N
  ),
  Gene
]


cat(
  "\n============================================================\n",
  "RA ALE GENES\n",
  "============================================================\n",
  "LightGBM:\n",
  paste(
    LGBM_PROFILE_GENES,
    collapse = "\n"
  ),
  "\n\nTabPFN:\n",
  paste(
    TABPFN_PROFILE_GENES,
    collapse = "\n"
  ),
  "\n============================================================\n",
  sep = ""
)


ale_gene_table <- data.table::rbindlist(
  list(

    data.table::data.table(
      Model = "LightGBM",
      Rank = seq_along(
        LGBM_PROFILE_GENES
      ),
      Gene = LGBM_PROFILE_GENES
    ),

    data.table::data.table(
      Model = "TabPFN",
      Rank = seq_along(
        TABPFN_PROFILE_GENES
      ),
      Gene = TABPFN_PROFILE_GENES
    )
  )
)


data.table::fwrite(
  ale_gene_table,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_ALE_selected_genes.csv"
  )
)


# ======================================================================
# 15. CONFIRM LIVE EXPLAINERS ARE STILL PRESENT
#
# Your current R session should contain these from the previous run.
# ======================================================================

required_ale_objects <- c(
  "lgbm_explainer",
  "tabpfn_explainer"
)


missing_ale_objects <- required_ale_objects[
  !vapply(
    required_ale_objects,
    exists,
    inherits = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (
  length(
    missing_ale_objects
  ) > 0L
) {

  stop(
    "The current R session is missing:\n",
    paste(
      missing_ale_objects,
      collapse = "\n"
    ),
    "\n\nDo not rerun the models. Reload the saved final models and recreate ",
    "the explainers on RA_foundation_final_x_ice.rds."
  )
}


# ======================================================================
# 16. SUPPRESS TABPFN PYTORCH WARNING NOISE
# ======================================================================

reticulate::py_run_string("
import warnings
warnings.simplefilter('ignore', UserWarning)
")


# ======================================================================
# 17. LIGHTGBM ALE
# ======================================================================

cat(
  "\n============================================================\n",
  "CALCULATING LIGHTGBM ALE\n",
  "============================================================\n"
)


set.seed(
  7832L
)


lgbm_ale <- ingredients::accumulated_dependency(
  lgbm_explainer,
  variables = LGBM_PROFILE_GENES,
  N = NULL
)


message(
  "PASS: LightGBM ALE complete."
)


# ======================================================================
# 18. SAVE LIGHTGBM ALE IMMEDIATELY
# ======================================================================

saveRDS(
  lgbm_ale,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_ALE.rds"
  )
)


data.table::fwrite(
  data.table::as.data.table(
    lgbm_ale
  ),
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_ALE.csv"
  )
)


# ======================================================================
# 19. TABPFN ALE
# ======================================================================

cat(
  "\n============================================================\n",
  "CALCULATING TABPFN ALE\n",
  "============================================================\n"
)


set.seed(
  7832L
)


tabpfn_ale <- ingredients::accumulated_dependency(
  tabpfn_explainer,
  variables = TABPFN_PROFILE_GENES,
  N = NULL
)


message(
  "PASS: TabPFN ALE complete."
)


# ======================================================================
# 20. SAVE TABPFN ALE IMMEDIATELY
# ======================================================================

saveRDS(
  tabpfn_ale,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_TabPFN_ALE.rds"
  )
)


data.table::fwrite(
  data.table::as.data.table(
    tabpfn_ale
  ),
  file.path(
    DALEX_COMPARE_DIR,
    "RA_TabPFN_ALE.csv"
  )
)


# ======================================================================
# 21. LIGHTGBM ALE PLOT
# ======================================================================

p_lgbm_ale <- plot(
  lgbm_ale
) +

  ggplot2::labs(
    title =
      "LightGBM expression-response profiles",
    subtitle =
      "Accumulated local effects for the six most influential RA genes",
    x =
      "Gene expression [log2(CPM + 1)]",
    y =
      "Accumulated local effect"
  ) +

  publication_theme() +

  ggplot2::theme(
    strip.text =
      ggplot2::element_text(
        face = "italic",
        colour = "black"
      )
  )


print(
  p_lgbm_ale
)


# ======================================================================
# 22. TABPFN ALE PLOT
# ======================================================================

p_tabpfn_ale <- plot(
  tabpfn_ale
) +

  ggplot2::labs(
    title =
      "TabPFN expression-response profiles",
    subtitle =
      "Accumulated local effects for the six most influential RA genes",
    x =
      "Gene expression [log2(CPM + 1)]",
    y =
      "Accumulated local effect"
  ) +

  publication_theme() +

  ggplot2::theme(
    strip.text =
      ggplot2::element_text(
        face = "italic",
        colour = "black"
      )
  )


print(
  p_tabpfn_ale
)


# ======================================================================
# 23. SAVE INDIVIDUAL ALE FIGURES
# ======================================================================

ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_LightGBM_ALE.pdf"
  ),
  plot = p_lgbm_ale,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_TabPFN_ALE.pdf"
  ),
  plot = p_tabpfn_ale,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_LightGBM_ALE.png"
  ),
  plot = p_lgbm_ale,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_TabPFN_ALE.png"
  ),
  plot = p_tabpfn_ale,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 24. COMBINED ALE FIGURE
# ======================================================================

p_ale_combined <- (

  p_lgbm_ale +

    p_tabpfn_ale

) +

  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) +

  patchwork::plot_annotation(

    title =
      "Expression-response profiles of the most influential RA genes",

    subtitle =
      paste0(
        "Accumulated local effects from LightGBM and TabPFN ",
        "in the identical diagnostic held-out population"
      ),

    tag_levels =
      "A"
  )


print(
  p_ale_combined
)


# ======================================================================
# 25. SAVE COMBINED ALE FIGURE
# ======================================================================

ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_ALE_LightGBM_vs_TabPFN.pdf"
  ),
  plot = p_ale_combined,
  width = 13.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_ALE_LightGBM_vs_TabPFN.png"
  ),
  plot = p_ale_combined,
  width = 13.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 26. FINAL RESULT TABLE
# ======================================================================

final_interpretability_summary <- data.table::data.table(

  Disease =
    "RA",

  LightGBM_baseline_AUROC =
    lgbm_importance_summary$Baseline_AUROC[[1L]],

  TabPFN_baseline_AUROC =
    tabpfn_importance_summary$Baseline_AUROC[[1L]],

  Importance_Spearman_rho =
    importance_spearman,

  Shared_top20_genes =
    length(
      shared_top20
    ),

  Top20_Jaccard =
    top20_jaccard,

  LightGBM_top_gene =
    lgbm_importance_summary$Gene[[1L]],

  TabPFN_top_gene =
    tabpfn_importance_summary$Gene[[1L]]
)


print(
  final_interpretability_summary
)


data.table::fwrite(
  final_interpretability_summary,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_interpretability_summary.csv"
  )
)


# ======================================================================
# 27. OUTPUT AUDIT
# ======================================================================

expected_files <- c(

  "RA_LightGBM_vs_TabPFN_importance_combined.csv",

  "Table_RA_top20_importance_LightGBM_vs_TabPFN.csv",

  "Table_RA_importance_concordance.csv",

  "Table_RA_shared_top20_genes.csv",

  "Table_RA_ALE_selected_genes.csv",

  "Table_RA_interpretability_summary.csv",

  "RA_LightGBM_ALE.rds",

  "RA_TabPFN_ALE.rds",

  "RA_LightGBM_ALE.csv",

  "RA_TabPFN_ALE.csv",

  "Figure_RA_DALEX_LightGBM_vs_TabPFN.pdf",

  "Figure_RA_DALEX_LightGBM_vs_TabPFN.png",

  "Figure_RA_DALEX_importance_concordance.pdf",

  "Figure_RA_LightGBM_ALE.pdf",

  "Figure_RA_TabPFN_ALE.pdf",

  "Figure_RA_ALE_LightGBM_vs_TabPFN.pdf",

  "Figure_RA_ALE_LightGBM_vs_TabPFN.png"
)


output_audit <- data.table::data.table(

  File =
    expected_files,

  Exists =
    file.exists(
      file.path(
        DALEX_COMPARE_DIR,
        expected_files
      )
    )
)


print(
  output_audit
)


data.table::fwrite(
  output_audit,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_FINAL_DALEX_ALE_output_audit.csv"
  )
)


# ======================================================================
# 28. FINAL CONSOLE SUMMARY
# ======================================================================

cat(
  "\n\n============================================================\n",
  "RA LIGHTGBM vs TABPFN INTERPRETABILITY COMPLETE\n",
  "============================================================\n",
  "LightGBM baseline AUROC: ",
  sprintf(
    "%.3f",
    lgbm_importance_summary$Baseline_AUROC[[1L]]
  ),
  "\n",
  "TabPFN baseline AUROC: ",
  sprintf(
    "%.3f",
    tabpfn_importance_summary$Baseline_AUROC[[1L]]
  ),
  "\n",
  "Feature-importance Spearman rho: ",
  sprintf(
    "%.3f",
    importance_spearman
  ),
  "\n",
  "Shared top-20 genes: ",
  length(
    shared_top20
  ),
  "/20\n",
  "Top-20 Jaccard: ",
  sprintf(
    "%.3f",
    top20_jaccard
  ),
  "\n\n",
  "LightGBM ALE genes:\n",
  paste(
    LGBM_PROFILE_GENES,
    collapse = ", "
  ),
  "\n\n",
  "TabPFN ALE genes:\n",
  paste(
    TABPFN_PROFILE_GENES,
    collapse = ", "
  ),
  "\n\n",
  "Saved to:\n",
  DALEX_COMPARE_DIR,
  "\n============================================================\n",
  sep = ""
)
