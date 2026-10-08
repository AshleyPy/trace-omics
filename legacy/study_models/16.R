# ======================================================================
# PART 14
# DALEX GLOBAL VARIABLE IMPORTANCE
#
# OUTPUT:
#   - Full DALEX permutation results saved for all genes
#   - Publication figure shows top 20 genes only
#   - Mean permutation importance
#   - 95% empirical interval across permutations
#   - Legend positioned at top
#
# IMPORTANCE:
#   Increase in one-minus-AUROC loss after gene permutation
#
# Larger positive values indicate a greater contribution to
# discrimination.
# ======================================================================


# ======================================================================
# 14.1 LOSS FUNCTION
# ======================================================================

loss_one_minus_auc <- function(
  observed,
  predicted
) {

  truth_factor <- factor(
    ifelse(
      observed == 1,
      POSITIVE_CLASS,
      NEGATIVE_CLASS
    ),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  roc_object <- pROC::roc(
    response = truth_factor,
    predictor = predicted,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )

  1 -
    as.numeric(
      pROC::auc(
        roc_object
      )
    )
}


# ======================================================================
# 14.2 RUN DALEX PERMUTATION IMPORTANCE
# ======================================================================

DALEX_PERMUTATIONS <- 100L

TOP_DALEX_GENES <- 20L


set.seed(
  SEED_PUBLICATION
)


dalex_variable_importance <- DALEX::model_parts(
  tb_rf_explainer,
  loss_function = loss_one_minus_auc,
  type = "difference",
  B = DALEX_PERMUTATIONS
)


print(
  dalex_variable_importance
)


# ======================================================================
# 14.3 SAVE RAW DALEX OBJECT
# ======================================================================

saveRDS(
  dalex_variable_importance,
  file.path(
    DALEX_DIR,
    "TB_final_RF_DALEX_variable_importance.rds"
  )
)


dalex_variable_importance_dt <- data.table::as.data.table(
  dalex_variable_importance
)


data.table::fwrite(
  dalex_variable_importance_dt,
  file.path(
    DALEX_DIR,
    "TB_final_RF_DALEX_variable_importance_ALL.csv"
  )
)


# ======================================================================
# 14.4 REMOVE DALEX REFERENCE ROWS
# ======================================================================

dalex_gene_dt <- data.table::copy(
  dalex_variable_importance_dt
)


dalex_gene_dt <- dalex_gene_dt[
  !variable %in% c(
    "_baseline_",
    "_full_model_"
  )
]


if (
  nrow(
    dalex_gene_dt
  ) ==
    0L
) {

  stop(
    "No gene-level DALEX importance values were found."
  )
}


# ======================================================================
# 14.5 SUMMARISE IMPORTANCE ACROSS PERMUTATIONS
# ======================================================================

dalex_gene_summary <- dalex_gene_dt[
  ,
  .(

    Mean_importance =
      mean(
        dropout_loss,
        na.rm = TRUE
      ),

    Median_importance =
      stats::median(
        dropout_loss,
        na.rm = TRUE
      ),

    SD_importance =
      stats::sd(
        dropout_loss,
        na.rm = TRUE
      ),

    Lower95 =
      as.numeric(
        stats::quantile(
          dropout_loss,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

    Upper95 =
      as.numeric(
        stats::quantile(
          dropout_loss,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

    Permutations =
      sum(
        is.finite(
          dropout_loss
        )
      )
  ),
  by = .(
    variable,
    label
  )
]


# ======================================================================
# 14.6 RANK GENES
# ======================================================================

data.table::setorder(
  dalex_gene_summary,
  -Mean_importance,
  variable
)


dalex_gene_summary[
  ,
  Rank := seq_len(
    .N
  )
]


data.table::fwrite(
  dalex_gene_summary,
  file.path(
    TABLE_DIR,
    "Table_S_DALEX_global_importance_all_genes.csv"
  )
)


print(
  dalex_gene_summary[
    1:min(
      25L,
      .N
    )
  ]
)


# ======================================================================
# 14.7 TOP GENES FOR PUBLICATION FIGURE
# ======================================================================

dalex_plot_dt <- data.table::copy(
  dalex_gene_summary[
    1:min(
      TOP_DALEX_GENES,
      .N
    )
  ]
)


# Most important gene at the top after coord_flip().
dalex_plot_dt[
  ,
  variable := factor(
    variable,
    levels = rev(
      variable
    )
  )
]


# Cleaner model label for legend.
dalex_plot_dt[
  ,
  Model := paste0(
    "logCPM + Random Forest | ",
    LOCKED_TOP_N,
    " genes"
  )
]


# ======================================================================
# 14.8 PUBLICATION FIGURE
# ======================================================================

p_dalex_global <- ggplot2::ggplot(
  dalex_plot_dt,
  ggplot2::aes(
    x = variable,
    y = Mean_importance,
    fill = Model
  )
) +

  # Zero = no deterioration in AUROC-related loss after permutation.
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.55
  ) +

  # Mean importance.
  ggplot2::geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.20
  ) +

  # 95% empirical permutation interval.
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = Lower95,
      ymax = Upper95
    ),
    width = 0.16,
    linewidth = 0.45
  ) +

  ggplot2::coord_flip() +

  ggplot2::scale_fill_manual(
    values = c(
      "logCPM + Random Forest | 100 genes" =
        SET1[[2L]],
      "logCPM + Random Forest | 50 genes" =
        SET1[[2L]],
      "logCPM + Random Forest | 25 genes" =
        SET1[[2L]]
    ),
    name = NULL
  ) +

  ggplot2::labs(
    title = "DALEX global feature importance",
    subtitle = paste0(
      "Top ",
      TOP_DALEX_GENES,
      " genes ranked by permutation-based increase in one-minus-AUROC loss"
    ),
    x = NULL,
    y = "Increase in one-minus-AUROC loss"
  ) +

  publication_theme() +

  ggplot2::theme(

    # --------------------------------------------------------------
    # Legend at top
    # --------------------------------------------------------------

    legend.position = "top",

    legend.justification = "left",

    legend.direction = "horizontal",

    # --------------------------------------------------------------
    # Gene labels
    # --------------------------------------------------------------

    axis.text.y = ggplot2::element_text(
      size = 9.5,
      face = "italic",
      colour = "black"
    ),

    # --------------------------------------------------------------
    # Improve readability
    # --------------------------------------------------------------

    panel.grid.major.y =
      ggplot2::element_blank(),

    panel.grid.minor =
      ggplot2::element_blank(),

    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 13
      ),

    plot.subtitle =
      ggplot2::element_text(
        size = 9.5
      ),

    plot.margin =
      ggplot2::margin(
        t = 8,
        r = 12,
        b = 8,
        l = 8
      )
  )


print(
  p_dalex_global
)


# ======================================================================
# 14.9 SAVE TOP-GENE PLOT DATA
# ======================================================================

data.table::fwrite(
  dalex_plot_dt,
  file.path(
    FIGURE_DIR,
    "Figure_S2_DALEX_global_importance_plot_data.csv"
  )
)


# ======================================================================
# 14.10 SAVE FIGURE
# ======================================================================

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_S2_DALEX_global_importance.pdf"
  ),
  plot = p_dalex_global,
  width = 8.2,
  height = 6.8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_S2_DALEX_global_importance.png"
  ),
  plot = p_dalex_global,
  width = 8.2,
  height = 6.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


cat(
  "\n============================================================\n",
  "DALEX GLOBAL IMPORTANCE COMPLETE\n",
  "============================================================\n",
  "Total genes evaluated: ",
  data.table::uniqueN(
    dalex_gene_dt$variable
  ),
  "\n",
  "Genes displayed: ",
  nrow(
    dalex_plot_dt
  ),
  "\n",
  "Permutations requested: ",
  DALEX_PERMUTATIONS,
  "\n",
  "Top gene: ",
  as.character(
    dalex_gene_summary$variable[[1L]]
  ),
  "\n",
  "Top mean importance: ",
  sprintf(
    "%.4f",
    dalex_gene_summary$Mean_importance[[1L]]
  ),
  "\n",
  "============================================================\n",
  sep = ""
)
