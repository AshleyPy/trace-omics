# ======================================================================
# STRICT LOCO
# LIGHTGBM vs TABPFN
#
# TARGET:
#   ACS
#   Catalysis
#   GC6
#
# PRIMARY:
#   AUROC
#
# MODEL SELECTION:
#   LightGBM inner CV on source cohorts only.
#
# TabPFN:
#   no task-specific hyperparameter tuning.
# ======================================================================

GENERALISATION_DIR <- file.path(
  FOUNDATION_DIR,
  "strict_LOCO_generalisation"
)

dir.create(
  GENERALISATION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


LOCO_TOP_N <- 100L

LOCO_INNER_FOLDS <- 3L

# ======================================================================
# SOURCE-ONLY LIGHTGBM TUNING
# ======================================================================

tune_lightgbm_source_only <- function(
  source_counts,
  source_metadata,
  target_index
) {

  inner_fold_id <- make_joint_stratified_folds(
    metadata = source_metadata,
    number_of_folds = LOCO_INNER_FOLDS,
    seed = SEED_FOUNDATION +
      target_index *
        100000L
  )


  source_metadata <- data.table::copy(
    source_metadata
  )


  source_metadata[
    ,
    loco_inner_fold := inner_fold_id
  ]


  tuning_rows <- list()

  tuning_counter <- 0L


  for (
    inner_fold in seq_len(
      LOCO_INNER_FOLDS
    )
  ) {

    inner_train_metadata <- source_metadata[
      loco_inner_fold !=
        inner_fold
    ]


    inner_valid_metadata <- source_metadata[
      loco_inner_fold ==
        inner_fold
    ]


    inner_train_counts <- source_counts[
      ,
      inner_train_metadata$sample_id,
      drop = FALSE
    ]


    inner_valid_counts <- source_counts[
      ,
      inner_valid_metadata$sample_id,
      drop = FALSE
    ]


    inner_representation <- fit_apply_logcpm(
      train_counts = inner_train_counts,
      new_counts = inner_valid_counts,
      train_metadata = inner_train_metadata
    )


    inner_genes <- select_rf_features_final(
      train_matrix = inner_representation$train,
      train_metadata = inner_train_metadata,
      top_n = LOCO_TOP_N,
      seed = SEED_FOUNDATION +
        target_index *
          1000000L +
        inner_fold *
          1000L
    )


    inner_x_train <- t(
      inner_representation$train[
        inner_genes,
        ,
        drop = FALSE
      ]
    )


    inner_x_valid <- t(
      inner_representation$new[
        inner_genes,
        ,
        drop = FALSE
      ]
    )


    for (
      config_index in seq_len(
        nrow(
          lightgbm_grid
        )
      )
    ) {

      current_config <- lightgbm_grid[
        config_index
      ]


      probability <- fit_predict_lightgbm(
        x_train = inner_x_train,
        y_train = inner_train_metadata$outcome,
        x_test = inner_x_valid,
        parameters = current_config,
        seed = SEED_FOUNDATION +
          target_index *
            10000000L +
          inner_fold *
            10000L +
          config_index
      )


      metric <- evaluate_binary_predictions(
        truth = inner_valid_metadata$outcome,
        probability = probability
      )


      tuning_counter <- tuning_counter +
        1L


      tuning_rows[[tuning_counter]] <- data.table::data.table(
        Inner_fold = inner_fold,
        Config_ID = current_config$Config_ID[[1L]],
        AUROC = metric$AUROC[[1L]]
      )
    }
  }


  tuning_dt <- data.table::rbindlist(
    tuning_rows
  )


  tuning_summary <- tuning_dt[
    ,
    .(
      Mean_inner_AUROC = mean(
        AUROC
      ),
      SD_inner_AUROC = stats::sd(
        AUROC
      )
    ),
    by = Config_ID
  ]


  data.table::setorder(
    tuning_summary,
    -Mean_inner_AUROC,
    SD_inner_AUROC,
    Config_ID
  )


  best_config_id <- tuning_summary$Config_ID[[1L]]


  list(
    best_config =
      lightgbm_grid[
        Config_ID ==
          best_config_id
      ],
    summary =
      tuning_summary,
    raw =
      tuning_dt
  )
}

# ======================================================================
# STRICT LOCO LOOP
# ======================================================================

LOCO_COHORT_ORDER <- c(
  "ACS",
  "Catalysis",
  "GC6"
)


loco_prediction_rows <- list()
loco_metric_rows <- list()
loco_tuning_rows <- list()
loco_gene_rows <- list()


prediction_counter <- 0L
metric_counter <- 0L
tuning_counter <- 0L
gene_counter <- 0L


for (
  target_index in seq_along(
    LOCO_COHORT_ORDER
  )
) {

  target_cohort <- LOCO_COHORT_ORDER[[target_index]]


  cat(
    "\n\n============================================================\n",
    "STRICT LOCO: ",
    target_cohort,
    " COMPLETELY HELD OUT\n",
    "============================================================\n",
    sep = ""
  )


  source_metadata <- development_metadata[
    cohort_group !=
      target_cohort
  ]


  target_metadata <- development_metadata[
    cohort_group ==
      target_cohort
  ]


  source_counts <- development_counts[
    ,
    source_metadata$sample_id,
    drop = FALSE
  ]


  target_counts <- development_counts[
    ,
    target_metadata$sample_id,
    drop = FALSE
  ]


  # ====================================================================
  # TUNE LIGHTGBM ON SOURCE COHORTS ONLY
  # ====================================================================

  lightgbm_tuning <- tune_lightgbm_source_only(
    source_counts = source_counts,
    source_metadata = source_metadata,
    target_index = target_index
  )


  best_lgbm_config <- lightgbm_tuning$best_config


  cat(
    "Selected LightGBM config: ",
    best_lgbm_config$Config_ID[[1L]],
    "\n",
    sep = ""
  )


  tuning_counter <- tuning_counter +
    1L


  loco_tuning_rows[[tuning_counter]] <- data.table::copy(
    lightgbm_tuning$summary
  )[
    ,
    Target_cohort := target_cohort
  ]


  # ====================================================================
  # REPRESENTATION
  # ====================================================================

  loco_representation <- fit_apply_logcpm(
    train_counts = source_counts,
    new_counts = target_counts,
    train_metadata = source_metadata
  )


  # ====================================================================
  # SAME FEATURES FOR BOTH MODELS
  # ====================================================================

  loco_genes <- select_rf_features_final(
    train_matrix = loco_representation$train,
    train_metadata = source_metadata,
    top_n = LOCO_TOP_N,
    seed = SEED_FOUNDATION +
      target_index *
        20000000L
  )


  loco_x_train <- t(
    loco_representation$train[
      loco_genes,
      ,
      drop = FALSE
    ]
  )


  loco_x_test <- t(
    loco_representation$new[
      loco_genes,
      ,
      drop = FALSE
    ]
  )


  gene_counter <- gene_counter +
    1L


  loco_gene_rows[[gene_counter]] <- data.table::data.table(
    Target_cohort = target_cohort,
    Rank = seq_along(
      loco_genes
    ),
    Gene = loco_genes
  )


  # ====================================================================
  # LIGHTGBM
  # ====================================================================

  lgbm_probability <- fit_predict_lightgbm(
    x_train = loco_x_train,
    y_train = source_metadata$outcome,
    x_test = loco_x_test,
    parameters = best_lgbm_config,
    seed = SEED_FOUNDATION +
      target_index *
        30000000L
  )


  # ====================================================================
  # TABPFN
  # ====================================================================

  tabpfn_probability <- fit_predict_tabpfn(
    x_train = loco_x_train,
    y_train = source_metadata$outcome,
    x_test = loco_x_test,
    seed = SEED_FOUNDATION +
      target_index *
        40000000L
  )


  # ====================================================================
  # EVALUATE BOTH
  # ====================================================================

  for (
    current_model in MODEL_ORDER
  ) {

    probability <- if (
      current_model ==
        "LightGBM"
    ) {
      lgbm_probability
    } else {
      tabpfn_probability
    }


    current_roc <- pROC::roc(
      response = factor(
        target_metadata$outcome,
        levels = c(
          NEGATIVE_CLASS,
          POSITIVE_CLASS
        )
      ),
      predictor = probability,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      ),
      direction = "<",
      quiet = TRUE
    )


    current_ci <- as.numeric(
      pROC::ci.auc(
        current_roc,
        method = "delong",
        conf.level = 0.95
      )
    )


    current_metrics <- evaluate_binary_predictions(
      truth = target_metadata$outcome,
      probability = probability
    )


    prediction_counter <- prediction_counter +
      1L


    loco_prediction_rows[[prediction_counter]] <- data.table::data.table(
      Target_cohort = target_cohort,
      Model = current_model,
      sample_id = target_metadata$sample_id,
      truth = as.character(
        target_metadata$outcome
      ),
      probability_case = probability
    )


    metric_counter <- metric_counter +
      1L


    loco_metric_rows[[metric_counter]] <- data.table::data.table(
      Target_cohort = target_cohort,
      Model = current_model,
      n = nrow(
        target_metadata
      ),
      AUROC = as.numeric(
        pROC::auc(
          current_roc
        )
      ),
      AUROC_lower95 = current_ci[[1L]],
      AUROC_upper95 = current_ci[[3L]],
      AUPRC = current_metrics$AUPRC[[1L]],
      BACC = current_metrics$BACC[[1L]],
      Sensitivity = current_metrics$Sensitivity[[1L]],
      Specificity = current_metrics$Specificity[[1L]],
      Brier = current_metrics$Brier[[1L]]
    )
  }
}

# ======================================================================
# GENERALISABILITY SUMMARY
# ======================================================================

loco_foundation_predictions <- rbindlist(
  loco_prediction_rows
)


loco_foundation_metrics <- rbindlist(
  loco_metric_rows
)


loco_foundation_tuning <- rbindlist(
  loco_tuning_rows
)


loco_foundation_genes <- rbindlist(
  loco_gene_rows
)


generalisability_summary <- loco_foundation_metrics[
  ,
  .(
    Mean_LOCO_AUROC = mean(
      AUROC
    ),
    Worst_cohort_AUROC = min(
      AUROC
    ),
    SD_LOCO_AUROC = sd(
      AUROC
    ),
    Mean_AUPRC = mean(
      AUPRC
    ),
    Mean_BACC = mean(
      BACC
    ),
    Mean_Brier = mean(
      Brier
    )
  ),
  by = Model
]


data.table::setorder(
  generalisability_summary,
  -Worst_cohort_AUROC,
  -Mean_LOCO_AUROC
)


generalisability_summary[
  ,
  Generalisability_rank :=
    seq_len(
      .N
    )
]


print(
  loco_foundation_metrics
)


print(
  generalisability_summary
)

# ======================================================================
# GENERALISATION / TRANSPORTABILITY SHIFT
# ======================================================================

nested_reference <- foundation_pooled_summary[
  ,
  .(
    Model,
    Nested_OOF_AUROC = AUROC
  )
]


transportability_table <- merge(
  generalisability_summary,
  nested_reference,
  by = "Model"
)


transportability_table[
  ,
  Delta_LOCO_minus_nested :=
    Mean_LOCO_AUROC -
      Nested_OOF_AUROC
]


transportability_table[
  ,
  AUROC_retention_percent :=
    100 *
      Mean_LOCO_AUROC /
      Nested_OOF_AUROC
]


data.table::setorder(
  transportability_table,
  -Worst_cohort_AUROC
)


print(
  transportability_table
)

# ======================================================================
# COHORT-SPECIFIC GENERALISABILITY FIGURE
# ======================================================================

loco_plot_dt <- data.table::copy(
  loco_foundation_metrics
)


loco_plot_dt[
  ,
  Target_cohort := factor(
    Target_cohort,
    levels = LOCO_COHORT_ORDER
  )
]


loco_plot_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]


p_generalisation <- ggplot2::ggplot(
  loco_plot_dt,
  ggplot2::aes(
    x = Target_cohort,
    y = AUROC,
    colour = Model,
    group = Model
  )
) +

  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey45"
  ) +

  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = AUROC_lower95,
      ymax = AUROC_upper95
    ),
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    width = 0.12,
    linewidth = 0.60
  ) +

  ggplot2::geom_line(
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    linewidth = 0.80
  ) +

  ggplot2::geom_point(
    position = ggplot2::position_dodge(
      width = 0.35
    ),
    size = 3.2
  ) +

  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    name = "Model"
  ) +

  ggplot2::coord_cartesian(
    ylim = c(
      0.45,
      1
    )
  ) +

  ggplot2::labs(
    title = "Transportability to completely unseen TB cohorts",
    subtitle = "Each point is obtained after excluding the entire target cohort from preprocessing, feature selection and model fitting",
    x = "Completely held-out cohort",
    y = "AUROC (95% CI)"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "top"
  )


print(
  p_generalisation
)


ggplot2::ggsave(
  filename = file.path(
    GENERALISATION_DIR,
    "Figure_STRICT_LOCO_LightGBM_vs_TabPFN.pdf"
  ),
  plot = p_generalisation,
  width = 7.8,
  height = 5.5,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)
data.table::fwrite(
  loco_foundation_predictions,
  file.path(
    GENERALISATION_DIR,
    "01_strict_LOCO_predictions.csv"
  )
)


data.table::fwrite(
  loco_foundation_metrics,
  file.path(
    GENERALISATION_DIR,
    "02_strict_LOCO_by_cohort_metrics.csv"
  )
)


data.table::fwrite(
  generalisability_summary,
  file.path(
    GENERALISATION_DIR,
    "03_generalisability_summary.csv"
  )
)


data.table::fwrite(
  transportability_table,
  file.path(
    GENERALISATION_DIR,
    "04_nestedCV_to_LOCO_retention.csv"
  )
)


data.table::fwrite(
  loco_foundation_tuning,
  file.path(
    GENERALISATION_DIR,
    "05_source_only_LightGBM_tuning.csv"
  )
)


data.table::fwrite(
  loco_foundation_genes,
  file.path(
    GENERALISATION_DIR,
    "06_LOCO_selected_genes.csv"
  )
)
