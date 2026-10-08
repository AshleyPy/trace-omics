# ======================================================================
# PART 9A
# DIAGNOSTIC ICE ROC
# EMPIRICAL STEPWISE ROC + 95% CI RIBBON
# ======================================================================

ice_truth_factor <- factor(
  diagnostic_ice_predictions$truth,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  )
)


ice_roc <- pROC::roc(
  response = ice_truth_factor,
  predictor = diagnostic_ice_predictions$probability_case,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


ice_auc <- as.numeric(
  pROC::auc(
    ice_roc
  )
)


ice_auc_ci <- as.numeric(
  pROC::ci.auc(
    ice_roc,
    method = "delong",
    conf.level = 0.95
  )
)


# ======================================================================
# EMPIRICAL COORDINATES
# ======================================================================

ice_coordinates <- pROC::coords(
  ice_roc,
  x = "all",
  ret = c(
    "specificity",
    "sensitivity"
  ),
  transpose = FALSE
)


ice_roc_dt <- data.table::data.table(

  False_positive_rate =
    1 -
    as.numeric(
      ice_coordinates$specificity
    ),

  True_positive_rate =
    as.numeric(
      ice_coordinates$sensitivity
    )
)


data.table::setorder(
  ice_roc_dt,
  False_positive_rate,
  True_positive_rate
)


# ======================================================================
# POINTWISE BOOTSTRAP CI
# ======================================================================

set.seed(
  SEED_PUBLICATION
)


ice_curve_ci <- pROC::ci.se(
  ice_roc,
  specificities = ROC_SPECIFICITY_GRID,
  conf.level = 0.95,
  boot.n = ROC_BOOTSTRAPS,
  boot.stratified = TRUE
)


ice_ribbon_dt <- data.table::data.table(

  False_positive_rate =
    1 -
    ROC_SPECIFICITY_GRID,

  Sensitivity_lower =
    as.numeric(
      ice_curve_ci[
        ,
        1L
      ]
    ),

  Sensitivity_upper =
    as.numeric(
      ice_curve_ci[
        ,
        3L
      ]
    )
)


data.table::setorder(
  ice_ribbon_dt,
  False_positive_rate
)


ice_roc_label <- sprintf(
  "AUROC %.3f (%.3f\u2013%.3f)",
  ice_auc,
  ice_auc_ci[[1L]],
  ice_auc_ci[[3L]]
)


# ======================================================================
# PLOT
# ======================================================================

p_ice_roc <- ggplot2::ggplot() +

  ggplot2::geom_ribbon(
    data = ice_ribbon_dt,
    ggplot2::aes(
      x = False_positive_rate,
      ymin = Sensitivity_lower,
      ymax = Sensitivity_upper
    ),
    fill = SET1[[1L]],
    alpha = 0.15
  ) +

  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.70
  ) +

  ggplot2::geom_step(
    data = ice_roc_dt,
    ggplot2::aes(
      x = False_positive_rate,
      y = True_positive_rate
    ),
    colour = SET1[[1L]],
    linewidth = 1.15,
    direction = "vh"
  ) +

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
    expand = c(
      0,
      0
    )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(
    title = "Diagnostic held-out test ROC",
    subtitle = ice_roc_label,
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none"
  )
