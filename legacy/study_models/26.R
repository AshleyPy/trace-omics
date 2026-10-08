# ======================================================================
# LIGHTGBM vs TABPFN
# PUBLICATION FIGURES FROM SHARED NESTED-CV OOF PREDICTIONS
# ======================================================================

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

library(data.table)
library(ggplot2)
library(pROC)
library(PRROC)
library(patchwork)


FOUNDATION_LABEL <- if (
  exists(
    "TABPFN_LABEL"
  )
) {
  TABPFN_LABEL
} else {
  "TabPFN"
}


MODEL_ORDER <- c(
  "LightGBM",
  FOUNDATION_LABEL
)


MODEL_COLOURS <- stats::setNames(
  SET1[
    1:2
  ],
  MODEL_ORDER
)


COMPARISON_DIR <- file.path(
  FOUNDATION_DIR,
  "publication"
)

dir.create(
  COMPARISON_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# CHUNK 1 — FIXED SETUP
# LIGHTGBM vs TABPFN PUBLICATION COMPARISON
#
# IMPORTANT:
# Reload the saved nested-CV predictions so that any factor conversion
# performed interactively cannot have destroyed the original model label.
# ======================================================================


# ======================================================================
# 1. RELOAD PRISTINE NESTED-CV RESULTS
# ======================================================================

foundation_predictions <- data.table::fread(
  file.path(
    FOUNDATION_DIR,
    "01_nestedCV_predictions.csv"
  )
)

foundation_fold_metrics <- data.table::fread(
  file.path(
    FOUNDATION_DIR,
    "02_nestedCV_fold_metrics.csv"
  )
)


# ======================================================================
# 2. INSPECT ORIGINAL MODEL LABELS
# ======================================================================

cat(
  "\nOriginal prediction model labels:\n"
)

print(
  unique(
    foundation_predictions$Model
  )
)

cat(
  "\nOriginal fold-metric model labels:\n"
)

print(
  unique(
    foundation_fold_metrics$Model
  )
)


# ======================================================================
# 3. CANONICAL PUBLICATION LABEL
#
# Your installed TabPFN is version 0.1.11, so call it "TabPFN"
# rather than "TabPFN-3".
# ======================================================================

FOUNDATION_LABEL <- "TabPFN"

MODEL_ORDER <- c(
  "LightGBM",
  FOUNDATION_LABEL
)


# ======================================================================
# 4. STANDARDISE MODEL NAMES
#
# This also catches earlier strings such as:
#   TabPFN-3
#   TabPFN
# ======================================================================

normalise_foundation_model_name <- function(
  x
) {

  x <- as.character(
    x
  )

  x[
    grepl(
      "^TabPFN",
      x,
      ignore.case = TRUE
    )
  ] <- FOUNDATION_LABEL

  x
}


foundation_predictions[
  ,
  Model := normalise_foundation_model_name(
    Model
  )
]


foundation_fold_metrics[
  ,
  Model := normalise_foundation_model_name(
    Model
  )
]


# ======================================================================
# 5. HARD AUDIT BEFORE FACTORISATION
# ======================================================================

prediction_audit <- foundation_predictions[
  ,
  .(
    Rows = .N,
    Unique_samples = data.table::uniqueN(
      sample_id
    )
  ),
  by = Model
]


print(
  prediction_audit
)


expected_development_n <- nrow(
  development_metadata
)


if (
  !setequal(
    prediction_audit$Model,
    MODEL_ORDER
  )
) {

  stop(
    paste0(
      "Unexpected model labels remain.\n",
      "Expected: ",
      paste(
        MODEL_ORDER,
        collapse = ", "
      ),
      "\nObserved: ",
      paste(
        prediction_audit$Model,
        collapse = ", "
      )
    )
  )
}


if (
  any(
    prediction_audit$Rows !=
      expected_development_n
  )
) {

  stop(
    "Each model should contain ",
    expected_development_n,
    " OOF predictions."
  )
}


if (
  any(
    prediction_audit$Unique_samples !=
      expected_development_n
  )
) {

  stop(
    "Each model should contain one prediction per development sample."
  )
}


if (
  anyNA(
    foundation_predictions$Model
  )
) {

  stop(
    "NA model labels remain after repair."
  )
}


message(
  "PASS: LightGBM and TabPFN each contain one OOF prediction per development sample."
)


lgbm_oof <- foundation_predictions[
  Model == "LightGBM"
][
  order(
    sample_id
  )
]


tabpfn_oof <- foundation_predictions[
  Model == FOUNDATION_LABEL
][
  order(
    sample_id
  )
]


stopifnot(
  identical(
    lgbm_oof$sample_id,
    tabpfn_oof$sample_id
  )
)


stopifnot(
  identical(
    lgbm_oof$truth,
    tabpfn_oof$truth
  )
)

# ======================================================================
# COMPARISON ROC
# ======================================================================

MODEL_ROC_BOOTSTRAPS <- 2000L

ROC_SPECIFICITY_GRID_COMPARE <- seq(
  0,
  1,
  length.out = 101L
)


roc_curve_rows <- list()
roc_ribbon_rows <- list()
roc_summary_rows <- list()
roc_labels <- character()


for (
  model_index in seq_along(
    MODEL_ORDER
  )
) {

  current_model <- MODEL_ORDER[[model_index]]

  current_dt <- foundation_predictions[
    Model == current_model
  ]

  current_roc <- pROC::roc(
    response = factor(
      current_dt$truth,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      )
    ),
    predictor = current_dt$probability_case,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )


  current_auc <- as.numeric(
    pROC::auc(
      current_roc
    )
  )


  current_ci <- as.numeric(
    pROC::ci.auc(
      current_roc,
      method = "delong",
      conf.level = 0.95
    )
  )


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
    Model = current_model,
    FPR = 1 - as.numeric(
      current_coordinates$specificity
    ),
    TPR = as.numeric(
      current_coordinates$sensitivity
    )
  )


  data.table::setorder(
    current_curve_dt,
    FPR,
    TPR
  )


  roc_curve_rows[[model_index]] <- current_curve_dt


  set.seed(
    SEED_FOUNDATION +
      model_index *
        10000L
  )


  current_curve_ci <- pROC::ci.se(
    current_roc,
    specificities = ROC_SPECIFICITY_GRID_COMPARE,
    conf.level = 0.95,
    boot.n = MODEL_ROC_BOOTSTRAPS,
    boot.stratified = TRUE
  )


  current_ribbon_dt <- data.table::data.table(
    Model = current_model,
    FPR = 1 - ROC_SPECIFICITY_GRID_COMPARE,
    Lower95 = as.numeric(
      current_curve_ci[
        ,
        1L
      ]
    ),
    Upper95 = as.numeric(
      current_curve_ci[
        ,
        3L
      ]
    )
  )


  data.table::setorder(
    current_ribbon_dt,
    FPR
  )


  roc_ribbon_rows[[model_index]] <- current_ribbon_dt


  roc_labels[[current_model]] <- paste0(
    current_model,
    ": ",
    sprintf(
      "%.3f",
      current_auc
    ),
    " (",
    sprintf(
      "%.3f",
      current_ci[[1L]]
    ),
    "\u2013",
    sprintf(
      "%.3f",
      current_ci[[3L]]
    ),
    ")"
  )


  roc_summary_rows[[model_index]] <- data.table::data.table(
    Model = current_model,
    AUROC = current_auc,
    AUROC_lower95 = current_ci[[1L]],
    AUROC_upper95 = current_ci[[3L]]
  )
}


comparison_roc_curve_dt <- rbindlist(
  roc_curve_rows
)

comparison_roc_ribbon_dt <- rbindlist(
  roc_ribbon_rows
)

comparison_roc_summary <- rbindlist(
  roc_summary_rows
)


comparison_roc_curve_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

comparison_roc_ribbon_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

# ======================================================================
# PAIRED DELONG TEST
# ======================================================================

lgbm_roc <- pROC::roc(
  response = factor(
    lgbm_oof$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),
  predictor = lgbm_oof$probability_case,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


tabpfn_roc <- pROC::roc(
  response = factor(
    tabpfn_oof$truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  ),
  predictor = tabpfn_oof$probability_case,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


paired_delong <- pROC::roc.test(
  lgbm_roc,
  tabpfn_roc,
  method = "delong",
  paired = TRUE
)


lgbm_auc <- as.numeric(
  pROC::auc(
    lgbm_roc
  )
)

tabpfn_auc <- as.numeric(
  pROC::auc(
    tabpfn_roc
  )
)


delta_auc <- (
  tabpfn_auc -
    lgbm_auc
)


delong_p <- as.numeric(
  paired_delong$p.value
)


comparison_annotation <- paste0(
  "\u0394AUROC (",
  FOUNDATION_LABEL,
  " - LightGBM) = ",
  sprintf(
    "%+.3f",
    delta_auc
  ),
  "\n",
  ifelse(
    delong_p < 0.001,
    "Paired DeLong p < 0.001",
    paste0(
      "Paired DeLong p = ",
      sprintf(
        "%.3f",
        delong_p
      )
    )
  )
)

p_model_roc <- ggplot2::ggplot() +

  ggplot2::geom_ribbon(
    data = comparison_roc_ribbon_dt,
    ggplot2::aes(
      x = FPR,
      ymin = Lower95,
      ymax = Upper95,
      fill = Model,
      group = Model
    ),
    alpha = 0.13,
    colour = NA
  ) +

  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.65
  ) +

  ggplot2::geom_step(
    data = comparison_roc_curve_dt,
    ggplot2::aes(
      x = FPR,
      y = TPR,
      colour = Model,
      group = Model
    ),
    direction = "vh",
    linewidth = 1.15
  ) +

  ggplot2::annotate(
    geom = "label",
    x = 0.97,
    y = 0.05,
    label = comparison_annotation,
    hjust = 1,
    vjust = 0,
    size = 3.2,
    fill = "white"
  ) +

  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    labels = roc_labels[
      MODEL_ORDER
    ],
    name = "Model \u2014 AUROC (95% CI)"
  ) +

  ggplot2::scale_fill_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    guide = "none"
  ) +

  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
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
      by = 0.2
    ),
    expand = c(
      0,
      0
    )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(
    title = "Nested-CV discrimination",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "top",
    legend.justification = "left"
  )


print(
  p_model_roc
)

# ======================================================================
# PR COMPARISON
# ======================================================================

pr_curve_rows <- list()
pr_ribbon_rows <- list()
pr_labels <- character()


comparison_prevalence <- mean(
  lgbm_oof$truth ==
    POSITIVE_CLASS
)


for (
  model_index in seq_along(
    MODEL_ORDER
  )
) {

  current_model <- MODEL_ORDER[[model_index]]

  current_dt <- foundation_predictions[
    Model == current_model
  ]

  current_pr <- PRROC::pr.curve(
    scores.class0 =
      current_dt$probability_case[
        current_dt$truth ==
          POSITIVE_CLASS
      ],
    scores.class1 =
      current_dt$probability_case[
        current_dt$truth ==
          NEGATIVE_CLASS
      ],
    curve = TRUE
  )


  pr_curve_rows[[model_index]] <- data.table::data.table(
    Model = current_model,
    Recall = as.numeric(
      current_pr$curve[
        ,
        1L
      ]
    ),
    Precision = as.numeric(
      current_pr$curve[
        ,
        2L
      ]
    )
  )


  current_ribbon <- bootstrap_pr_ci(
    truth = current_dt$truth,
    probability = current_dt$probability_case,
    recall_grid = PR_RECALL_GRID,
    bootstraps = 2000L,
    seed = SEED_FOUNDATION +
      model_index *
        20000L
  )


  current_ribbon[
    ,
    Model := current_model
  ]


  pr_ribbon_rows[[model_index]] <- current_ribbon


  pr_labels[[current_model]] <- paste0(
    current_model,
    ": AUPRC ",
    sprintf(
      "%.3f",
      current_pr$auc.integral
    )
  )
}


comparison_pr_dt <- rbindlist(
  pr_curve_rows
)

comparison_pr_ribbon_dt <- rbindlist(
  pr_ribbon_rows
)


comparison_pr_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]

comparison_pr_ribbon_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]


p_model_pr <- ggplot2::ggplot() +

  ggplot2::geom_ribbon(
    data = comparison_pr_ribbon_dt,
    ggplot2::aes(
      x = Recall,
      ymin = Precision_lower,
      ymax = Precision_upper,
      fill = Model,
      group = Model
    ),
    alpha = 0.13,
    colour = NA
  ) +

  ggplot2::geom_hline(
    yintercept = comparison_prevalence,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +

  ggplot2::geom_path(
    data = comparison_pr_dt,
    ggplot2::aes(
      x = Recall,
      y = Precision,
      colour = Model,
      group = Model
    ),
    linewidth = 1.10
  ) +

  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    labels = pr_labels[
      MODEL_ORDER
    ],
    name = NULL
  ) +

  ggplot2::scale_fill_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    guide = "none"
  ) +

  ggplot2::coord_cartesian(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  ggplot2::labs(
    title = "Nested-CV precision-recall",
    subtitle = paste0(
      "Dashed line = prevalence ",
      sprintf(
        "%.3f",
        comparison_prevalence
      )
    ),
    x = "Recall",
    y = "Precision"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "top"
  )

# ======================================================================
# CALIBRATION COMPARISON
# ======================================================================

comparison_calibration_rows <- list()


for (
  model_index in seq_along(
    MODEL_ORDER
  )
) {

  current_model <- MODEL_ORDER[[model_index]]

  current_dt <- foundation_predictions[
    Model == current_model
  ]


  comparison_calibration_rows[[model_index]] <- make_calibration_data(
    truth = current_dt$truth,
    probability = current_dt$probability_case,
    label = current_model,
    number_of_bins = 8L
  )
}


comparison_calibration_dt <- rbindlist(
  comparison_calibration_rows
)


comparison_calibration_dt[
  ,
  Group := factor(
    Group,
    levels = MODEL_ORDER
  )
]


data.table::setorder(
  comparison_calibration_dt,
  Group,
  Mean_predicted
)


p_model_calibration <- ggplot2::ggplot() +

  ggplot2::geom_ribbon(
    data = comparison_calibration_dt,
    ggplot2::aes(
      x = Mean_predicted,
      ymin = Lower95,
      ymax = Upper95,
      fill = Group,
      group = Group
    ),
    alpha = 0.13,
    colour = NA
  ) +

  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey40"
  ) +

  ggplot2::geom_line(
    data = comparison_calibration_dt,
    ggplot2::aes(
      x = Mean_predicted,
      y = Observed,
      colour = Group,
      group = Group
    ),
    linewidth = 0.95
  ) +

  ggplot2::geom_point(
    data = comparison_calibration_dt,
    ggplot2::aes(
      x = Mean_predicted,
      y = Observed,
      colour = Group
    ),
    size = 2.5
  ) +

  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    name = "Model"
  ) +

  ggplot2::scale_fill_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    guide = "none"
  ) +

  ggplot2::coord_equal(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  ggplot2::labs(
    title = "Nested-CV calibration",
    x = "Mean predicted probability",
    y = "Observed case proportion"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "top"
  )

# ======================================================================
# MAIN LIGHTGBM vs FOUNDATION-MODEL FIGURE
# ======================================================================

p_lgbm_vs_foundation <- (
  p_model_roc +
    p_model_pr +
    p_model_calibration
) +
  patchwork::plot_layout(
    ncol = 3
  ) +
  patchwork::plot_annotation(
    title = "Task-specific gradient boosting versus a tabular foundation model",
    subtitle = "Identical nested-CV splits, logCPM representation and training-selected transcriptomic features",
    tag_levels = "A"
  )


print(
  p_lgbm_vs_foundation
)


ggplot2::ggsave(
  filename = file.path(
    COMPARISON_DIR,
    "Figure_MAIN_LightGBM_vs_TabPFN_ROC_PR_calibration.pdf"
  ),
  plot = p_lgbm_vs_foundation,
  width = 15.8,
  height = 5.4,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    COMPARISON_DIR,
    "Figure_MAIN_LightGBM_vs_TabPFN_ROC_PR_calibration.png"
  ),
  plot = p_lgbm_vs_foundation,
  width = 15.8,
  height = 5.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ======================================================================
# PAIRED OUTER-FOLD AUROC FIGURE
# ======================================================================

paired_fold_dt <- data.table::copy(
  foundation_fold_metrics
)


paired_fold_dt[
  ,
  Model := factor(
    Model,
    levels = MODEL_ORDER
  )
]


p_paired_outer_auc <- ggplot2::ggplot(
  paired_fold_dt,
  ggplot2::aes(
    x = Model,
    y = AUROC,
    group = Outer_fold
  )
) +

  ggplot2::geom_line(
    colour = "grey65",
    linewidth = 0.70
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      colour = Model
    ),
    size = 3.2
  ) +

  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey45"
  ) +

  ggplot2::scale_colour_manual(
    values = MODEL_COLOURS,
    limits = MODEL_ORDER,
    guide = "none"
  ) +

  ggplot2::coord_cartesian(
    ylim = c(
      0.45,
      1
    )
  ) +

  ggplot2::labs(
    title = "Paired outer-fold discrimination",
    subtitle = "Lines connect predictions from the identical nested-CV outer fold",
    x = NULL,
    y = "Outer-fold AUROC"
  ) +

  publication_theme()


print(
  p_paired_outer_auc
)


ggplot2::ggsave(
  filename = file.path(
    COMPARISON_DIR,
    "Figure_S1_paired_outer_fold_AUROC.pdf"
  ),
  plot = p_paired_outer_auc,
  width = 6.4,
  height = 5.2,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)
