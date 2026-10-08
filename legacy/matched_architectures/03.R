# ======================================================================
# RA LIGHTGBM vs TABPFN
# FINAL DALEX IMPORTANCE FIGURE + ALE
#
# USES ONLY SAVED RESULTS FROM THE ALREADY-COMPLETED
# DIAGNOSTIC FINAL ICE ANALYSIS
#
# DOES NOT RERUN:
#   - preprocessing
#   - feature selection
#   - model tuning
#   - model fitting
#   - nested CV
#   - LOCO
#   - diagnostic ICE evaluation
#   - permutation importance
# ======================================================================


# ======================================================================
# 1. PATHS
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


ICE_DIR <- file.path(
  FOUNDATION_DIR,
  "diagnostic_FINAL_ICE"
)


MODEL_DIR <- file.path(
  FOUNDATION_DIR,
  "saved_models"
)


DALEX_COMPARE_DIR <- file.path(
  FOUNDATION_DIR,
  "DALEX_LightGBM_vs_TabPFN"
)


dir.create(
  DALEX_COMPARE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 2. R PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "RColorBrewer",
  "reticulate",
  "DALEX",
  "ingredients",
  "patchwork"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (
  length(
    missing_packages
  ) > 0L
) {

  stop(
    "Install these R packages first:\n",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


suppressPackageStartupMessages({

  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(reticulate)
  library(DALEX)
  library(ingredients)
  library(patchwork)

})


# ======================================================================
# 3. ENSURE THE CORRECT RA PYTHON ENVIRONMENT
#
# If the current R session already has the correct RA environment
# initialised, it will simply continue.
# ======================================================================

FOUNDATION_CONDA_ENV <- "ra_foundation_models_clean"


FOUNDATION_PYTHON <- reticulate::conda_python(
  FOUNDATION_CONDA_ENV
)


if (
  is.null(
    FOUNDATION_PYTHON
  ) ||
    !file.exists(
      FOUNDATION_PYTHON
    )
) {

  stop(
    "Could not resolve Python for:\n",
    FOUNDATION_CONDA_ENV
  )
}


requested_python <- normalizePath(
  FOUNDATION_PYTHON,
  winslash = "/",
  mustWork = TRUE
)


if (
  reticulate::py_available(
    initialize = FALSE
  )
) {

  active_python <- normalizePath(
    reticulate::py_config()$python,
    winslash = "/",
    mustWork = TRUE
  )


  if (
    !identical(
      tolower(
        active_python
      ),
      tolower(
        requested_python
      )
    )
  ) {

    stop(
      "The wrong Python environment is active.\n\n",
      "Active:\n",
      active_python,
      "\n\nRequired:\n",
      requested_python,
      "\n\nRestart R and rerun this script."
    )
  }

} else {

  Sys.setenv(
    RETICULATE_PYTHON =
      FOUNDATION_PYTHON
  )


  reticulate::use_python(
    FOUNDATION_PYTHON,
    required = TRUE
  )
}


np <- reticulate::import(
  "numpy",
  convert = FALSE
)


# Silence repeated PyTorch UserWarnings during ALE prediction.
reticulate::py_run_string("
import warnings
warnings.simplefilter('ignore', UserWarning)
")


cat(
  "\n============================================================\n",
  "RA PYTHON ENVIRONMENT READY\n",
  "============================================================\n",
  reticulate::py_config()$python,
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 4. LOAD SAVED RA EXPERIMENT STATE
# ======================================================================

state_file <- file.path(
  MODEL_DIR,
  "RA_LightGBM_TabPFN_experiment_state.rds"
)


if (
  !file.exists(
    state_file
  )
) {

  stop(
    "Could not find:\n",
    state_file
  )
}


foundation_state <- readRDS(
  state_file
)


POSITIVE_CLASS <- foundation_state$positive_class

NEGATIVE_CLASS <- foundation_state$negative_class

ICE_SELECTED_GENES <- foundation_state$selected_genes


stopifnot(
  identical(
    foundation_state$disease,
    "RA"
  )
)


cat(
  "\n============================================================\n",
  "RA FOUNDATION STATE LOADED\n",
  "============================================================\n",
  "Representation: ",
  foundation_state$representation,
  "\nSelected genes: ",
  length(
    ICE_SELECTED_GENES
  ),
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 5. LOAD THE EXACT SAVED DIAGNOSTIC ICE MATRIX
# ======================================================================

ice_matrix_file <- file.path(
  ICE_DIR,
  "RA_foundation_final_x_ice.rds"
)


if (
  !file.exists(
    ice_matrix_file
  )
) {

  stop(
    "Could not find:\n",
    ice_matrix_file
  )
}


ice_x_test <- readRDS(
  ice_matrix_file
)


ice_x_test <- as.matrix(
  ice_x_test
)


if (
  !identical(
    colnames(
      ice_x_test
    ),
    ICE_SELECTED_GENES
  )
) {

  stop(
    "Saved RA ICE matrix genes do not match the saved experiment state."
  )
}


# ======================================================================
# 6. LOAD ORIGINAL DIAGNOSTIC ICE PREDICTIONS
# ======================================================================

prediction_file <- file.path(
  ICE_DIR,
  "03_FINAL_ICE_predictions_LightGBM_vs_TabPFN.csv"
)


if (
  !file.exists(
    prediction_file
  )
) {

  stop(
    "Could not find:\n",
    prediction_file
  )
}


ice_predictions <- data.table::fread(
  prediction_file
)


ice_predictions[
  ,
  Model :=
    as.character(
      Model
    )
]


truth_dt <- unique(
  ice_predictions[
    ,
    .(
      sample_id,
      truth
    )
  ]
)


sample_match <- match(
  rownames(
    ice_x_test
  ),
  truth_dt$sample_id
)


if (
  anyNA(
    sample_match
  )
) {

  stop(
    "Could not align saved RA ICE samples to the original predictions."
  )
}


ice_truth <- truth_dt$truth[
  sample_match
]


dalex_compare_x <- as.data.frame(
  ice_x_test,
  check.names = FALSE
)


dalex_compare_y <- as.integer(
  ice_truth ==
    POSITIVE_CLASS
)


cat(
  "\n============================================================\n",
  "ORIGINAL RA DIAGNOSTIC ICE POPULATION\n",
  "============================================================\n",
  "Samples: ",
  nrow(
    dalex_compare_x
  ),
  "\nGenes: ",
  ncol(
    dalex_compare_x
  ),
  "\nCases: ",
  sum(
    ice_truth ==
      POSITIVE_CLASS
  ),
  "\nControls: ",
  sum(
    ice_truth ==
      NEGATIVE_CLASS
  ),
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 7. LOAD THE SAVED PERMUTATION-IMPORTANCE RESULTS
#
# IMPORTANT:
# These were already calculated.
# Nothing is recomputed here.
# ======================================================================

lgbm_importance_file <- file.path(
  DALEX_COMPARE_DIR,
  "RA_LightGBM_permutation_importance_summary.csv"
)


tabpfn_importance_file <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_summary.csv"
)


if (
  !file.exists(
    lgbm_importance_file
  )
) {

  stop(
    "Completed LightGBM importance summary was not found:\n",
    lgbm_importance_file
  )
}


if (
  !file.exists(
    tabpfn_importance_file
  )
) {

  stop(
    "Completed TabPFN importance summary was not found:\n",
    tabpfn_importance_file,
    "\n\nWait until the current TabPFN importance calculation finishes."
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
  Model :=
    "LightGBM"
]


tabpfn_importance_summary[
  ,
  Model :=
    "TabPFN"
]


data.table::setorder(
  lgbm_importance_summary,
  Rank
)


data.table::setorder(
  tabpfn_importance_summary,
  Rank
)


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
    "02_RA_LightGBM_TabPFN_permutation_importance_summary.csv"
  )
)


cat(
  "\n============================================================\n",
  "SAVED IMPORTANCE RESULTS LOADED\n",
  "============================================================\n",
  "LightGBM genes: ",
  nrow(
    lgbm_importance_summary
  ),
  "\nTabPFN genes: ",
  nrow(
    tabpfn_importance_summary
  ),
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
# 8. AESTHETICS
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
# 9. FUNCTION FOR THE TWO IMPORTANCE PANELS
#
# Same layout as the earlier LightGBM vs TabPFN comparison.
# ======================================================================

TOP_COMPARISON_GENES <- 20L


make_importance_plot <- function(
  importance_summary,
  model_name,
  top_n = 20L
) {

  plot_dt <- data.table::copy(
    importance_summary[
      Model ==
        model_name
    ][
      1:min(
        top_n,
        .N
      )
    ]
  )


  plot_dt[
    ,
    Gene :=
      factor(
        Gene,
        levels = rev(
          Gene
        )
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
      linewidth = 0.55
    ) +

    ggplot2::geom_col(
      width = 0.68,
      colour = "black",
      linewidth = 0.20
    ) +

    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Lower95_Delta_AUROC,
        ymax = Upper95_Delta_AUROC
      ),
      width = 0.16,
      linewidth = 0.45
    ) +

    ggplot2::coord_flip() +

    ggplot2::scale_fill_manual(
      values = MODEL_COLOURS,
      limits = MODEL_ORDER,
      name = NULL
    ) +

    ggplot2::labs(

      title =
        paste0(
          model_name,
          " global feature importance"
        ),

      subtitle =
        paste0(
          "Top ",
          top_n,
          " genes ranked by repeated permutation ",
          "\u0394AUROC"
        ),

      x =
        NULL,

      y =
        expression(
          "Permutation importance (" *
            Delta *
            "AUROC)"
        )
    ) +

    publication_theme() +

    ggplot2::theme(

      legend.position =
        "top",

      legend.justification =
        "left",

      axis.text.y =
        ggplot2::element_text(
          size = 9,
          face = "italic",
          colour = "black"
        ),

      panel.grid.major.y =
        ggplot2::element_blank()
    )
}


# ======================================================================
# 10. IMPORTANCE PANELS
# ======================================================================

p_lgbm_importance <- make_importance_plot(

  importance_summary =
    comparison_importance_summary,

  model_name =
    "LightGBM",

  top_n =
    TOP_COMPARISON_GENES
)


p_tabpfn_importance <- make_importance_plot(

  importance_summary =
    comparison_importance_summary,

  model_name =
    "TabPFN",

  top_n =
    TOP_COMPARISON_GENES
)


# ======================================================================
# 11. COMBINED DALEX IMPORTANCE FIGURE
# ======================================================================

p_model_importance_comparison <- (

  p_lgbm_importance +

    p_tabpfn_importance

) +

  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) +

  patchwork::plot_annotation(

    title =
      paste0(
        "Predictive feature reliance differs between ",
        "LightGBM and TabPFN in rheumatoid arthritis"
      ),

    subtitle =
      paste0(
        "Post-hoc permutation importance in the identical ",
        "diagnostic held-out population; 100 permutations per gene"
      ),

    tag_levels =
      "A"
  ) &

  ggplot2::theme(
    legend.position = "top"
  )


print(
  p_model_importance_comparison
)


# ======================================================================
# 12. SAVE IMPORTANCE FIGURE
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.pdf"
    ),

  plot =
    p_model_importance_comparison,

  width =
    13.5,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  device =
    grDevices::cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.png"
    ),

  plot =
    p_model_importance_comparison,

  width =
    13.5,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ======================================================================
# 13. SAVE TOP-20 IMPORTANCE TABLE
# ======================================================================

top20_importance_table <- data.table::rbindlist(
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
  top20_importance_table,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_top20_feature_importance_LightGBM_vs_TabPFN.csv"
  )
)


# ======================================================================
# 14. SELECT TOP SIX GENES FOR ALE
#
# Same method used previously:
# each model receives its own six highest-ranked genes.
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
  "GENES SELECTED FOR ALE\n",
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
# 15. LOAD SAVED FITTED MODELS
#
# Required only because ALE needs new predictions after locally changing
# gene-expression values.
#
# THESE MODELS ARE NOT REFITTED.
# ======================================================================

lgbm_model_file <- file.path(
  MODEL_DIR,
  "RA_FINAL_LightGBM_model.pickle"
)


tabpfn_model_file <- file.path(
  MODEL_DIR,
  "RA_FINAL_TabPFN_model.pickle"
)


if (
  !file.exists(
    lgbm_model_file
  )
) {

  stop(
    "Could not find:\n",
    lgbm_model_file
  )
}


if (
  !file.exists(
    tabpfn_model_file
  )
) {

  stop(
    "Could not find:\n",
    tabpfn_model_file
  )
}


FINAL_LGBM_MODEL <- reticulate::py_load_object(
  lgbm_model_file
)


FINAL_TABPFN_MODEL <- reticulate::py_load_object(
  tabpfn_model_file
)


# ======================================================================
# 16. PREDICTION FUNCTIONS FOR ALE
# ======================================================================

predict_lightgbm_dalex <- function(
  model,
  newdata
) {

  newdata_py <- np$array(
    as.matrix(
      newdata
    ),
    dtype = "float32"
  )


  probability_matrix <- reticulate::py_to_r(
    model$predict_proba(
      newdata_py
    )
  )


  as.numeric(
    probability_matrix[, 2L]
  )
}


predict_tabpfn_dalex <- function(
  model,
  newdata
) {

  newdata_py <- np$array(
    as.matrix(
      newdata
    ),
    dtype = "float32"
  )


  probability_matrix <- reticulate::py_to_r(
    model$predict_proba(
      newdata_py
    )
  )


  as.numeric(
    probability_matrix[, 2L]
  )
}


# ======================================================================
# 17. VERIFY MODELS STILL REPRODUCE ORIGINAL ICE PREDICTIONS
#
# This is an audit only.
# ======================================================================

audit_one_model <- function(
  model,
  model_name,
  predict_function
) {

  old_dt <- ice_predictions[
    Model ==
      model_name
  ]


  old_match <- match(
    rownames(
      ice_x_test
    ),
    old_dt$sample_id
  )


  old_probability <- old_dt$probability_case[
    old_match
  ]


  restored_probability <- predict_function(
    model,
    dalex_compare_x
  )


  data.table::data.table(

    Model =
      model_name,

    Max_absolute_difference =
      max(
        abs(
          old_probability -
            restored_probability
        )
      ),

    Mean_absolute_difference =
      mean(
        abs(
          old_probability -
            restored_probability
        )
      )
  )
}


prediction_audit <- data.table::rbindlist(
  list(

    audit_one_model(
      FINAL_LGBM_MODEL,
      "LightGBM",
      predict_lightgbm_dalex
    ),

    audit_one_model(
      FINAL_TABPFN_MODEL,
      "TabPFN",
      predict_tabpfn_dalex
    )
  )
)


print(
  prediction_audit
)


if (
  any(
    prediction_audit$Max_absolute_difference >
      1e-5
  )
) {

  stop(
    "The restored model does not reproduce the old diagnostic ICE ",
    "predictions. ALE has been stopped."
  )
}


data.table::fwrite(
  prediction_audit,
  file.path(
    DALEX_COMPARE_DIR,
    "00_RA_ALE_prediction_audit.csv"
  )
)


message(
  "PASS: saved models reproduce original RA diagnostic ICE predictions."
)


# ======================================================================
# 18. BUILD LIVE DALEX EXPLAINERS
#
# Same old ICE population.
# ======================================================================

lgbm_explainer <- DALEX::explain(

  model =
    FINAL_LGBM_MODEL,

  data =
    dalex_compare_x,

  y =
    dalex_compare_y,

  predict_function =
    predict_lightgbm_dalex,

  label =
    "LightGBM",

  type =
    "classification",

  verbose =
    FALSE
)


tabpfn_explainer <- DALEX::explain(

  model =
    FINAL_TABPFN_MODEL,

  data =
    dalex_compare_x,

  y =
    dalex_compare_y,

  predict_function =
    predict_tabpfn_dalex,

  label =
    "TabPFN",

  type =
    "classification",

  verbose =
    FALSE
)


# ======================================================================
# 19. LIGHTGBM ALE
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

  variables =
    LGBM_PROFILE_GENES,

  N =
    NULL
)


message(
  "PASS: LightGBM ALE complete."
)


# ======================================================================
# 20. TABPFN ALE
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

  variables =
    TABPFN_PROFILE_GENES,

  N =
    NULL
)


message(
  "PASS: TabPFN ALE complete."
)


# ======================================================================
# 21. SAVE ALE OBJECTS
# ======================================================================

saveRDS(
  lgbm_ale,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_ALE.rds"
  )
)


saveRDS(
  tabpfn_ale,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_TabPFN_ALE.rds"
  )
)


# ======================================================================
# 22. SAVE ALE DATA TABLES
# ======================================================================

lgbm_ale_dt <- data.table::as.data.table(
  lgbm_ale
)


tabpfn_ale_dt <- data.table::as.data.table(
  tabpfn_ale
)


data.table::fwrite(
  lgbm_ale_dt,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_ALE.csv"
  )
)


data.table::fwrite(
  tabpfn_ale_dt,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_TabPFN_ALE.csv"
  )
)


# ======================================================================
# 23. LIGHTGBM ALE FIGURE
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

    legend.position =
      "top",

    strip.text =
      ggplot2::element_text(
        face = "italic",
        size = 10,
        colour = "black"
      )
  )


print(
  p_lgbm_ale
)


# ======================================================================
# 24. TABPFN ALE FIGURE
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

    legend.position =
      "top",

    strip.text =
      ggplot2::element_text(
        face = "italic",
        size = 10,
        colour = "black"
      )
  )


print(
  p_tabpfn_ale
)


# ======================================================================
# 25. SAVE INDIVIDUAL ALE FIGURES
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_LightGBM_ALE.pdf"
    ),

  plot =
    p_lgbm_ale,

  width =
    9,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  device =
    grDevices::cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_TabPFN_ALE.pdf"
    ),

  plot =
    p_tabpfn_ale,

  width =
    9,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  device =
    grDevices::cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_LightGBM_ALE.png"
    ),

  plot =
    p_lgbm_ale,

  width =
    9,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_TabPFN_ALE.png"
    ),

  plot =
    p_tabpfn_ale,

  width =
    9,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ======================================================================
# 26. COMBINED ALE FIGURE
# ======================================================================

p_ale_comparison <- (

  p_lgbm_ale +

    p_tabpfn_ale

) +

  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) +

  patchwork::plot_annotation(

    title =
      "Expression-response profiles differ between LightGBM and TabPFN in rheumatoid arthritis",

    subtitle =
      paste0(
        "Accumulated local effects in the identical diagnostic ",
        "held-out population"
      ),

    tag_levels =
      "A"
  ) &

  ggplot2::theme(
    legend.position = "top"
  )


print(
  p_ale_comparison
)


# ======================================================================
# 27. SAVE COMBINED ALE FIGURE
# ======================================================================

ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_ALE_LightGBM_vs_TabPFN.pdf"
    ),

  plot =
    p_ale_comparison,

  width =
    13.5,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  device =
    grDevices::cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_ALE_LightGBM_vs_TabPFN.png"
    ),

  plot =
    p_ale_comparison,

  width =
    13.5,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ======================================================================
# 28. FEATURE-IMPORTANCE CONCORDANCE
#
# OPTIONAL BUT USEFUL FOR YOUR LightGBM vs TabPFN COMPARISON
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


importance_concordance <- data.table::data.table(

  Genes_compared =
    nrow(
      importance_wide
    ),

  Spearman_rho =
    importance_spearman,

  LightGBM_top20 =
    length(
      lgbm_top20
    ),

  TabPFN_top20 =
    length(
      tabpfn_top20
    ),

  Shared_top20 =
    length(
      shared_top20
    ),

  Top20_Jaccard =
    top20_jaccard
)


print(
  importance_concordance
)


data.table::fwrite(
  importance_concordance,
  file.path(
    DALEX_COMPARE_DIR,
    "Table_RA_feature_importance_concordance.csv"
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
# 29. CONCORDANCE FIGURE
# ======================================================================

importance_scatter_dt <- data.table::data.table(

  Gene =
    importance_wide$Gene,

  LightGBM =
    importance_wide$LightGBM,

  TabPFN =
    importance_wide$TabPFN
)


p_importance_concordance <- ggplot2::ggplot(

  importance_scatter_dt,

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
    size = 3.3
  ) +

  ggplot2::labs(

    title =
      "Concordance of RA transcriptomic feature reliance",

    subtitle =
      "Each point represents one gene evaluated in both models",

    x =
      expression(
        "LightGBM " *
          Delta *
          "AUROC"
      ),

    y =
      expression(
        "TabPFN " *
          Delta *
          "AUROC"
      )
  ) +

  publication_theme() +

  ggplot2::theme(
    legend.position = "none"
  )


print(
  p_importance_concordance
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_feature_importance_concordance.pdf"
    ),

  plot =
    p_importance_concordance,

  width =
    6.4,

  height =
    5.5,

  units =
    "in",

  dpi =
    600,

  device =
    grDevices::cairo_pdf
)


ggplot2::ggsave(

  filename =
    file.path(
      DALEX_COMPARE_DIR,
      "Figure_RA_feature_importance_concordance.png"
    ),

  plot =
    p_importance_concordance,

  width =
    6.4,

  height =
    5.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ======================================================================
# 30. FINAL OUTPUT AUDIT
# ======================================================================

expected_outputs <- c(

  "02_RA_LightGBM_TabPFN_permutation_importance_summary.csv",

  "Table_RA_top20_feature_importance_LightGBM_vs_TabPFN.csv",

  "Table_RA_ALE_selected_genes.csv",

  "RA_LightGBM_ALE.rds",

  "RA_TabPFN_ALE.rds",

  "RA_LightGBM_ALE.csv",

  "RA_TabPFN_ALE.csv",

  "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.pdf",

  "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.png",

  "Figure_RA_LightGBM_ALE.pdf",

  "Figure_RA_TabPFN_ALE.pdf",

  "Figure_RA_ALE_LightGBM_vs_TabPFN.pdf",

  "Figure_RA_ALE_LightGBM_vs_TabPFN.png",

  "Table_RA_feature_importance_concordance.csv",

  "Table_RA_shared_top20_genes.csv",

  "Figure_RA_feature_importance_concordance.pdf"
)


output_audit <- data.table::data.table(

  File =
    expected_outputs,

  Exists =
    file.exists(
      file.path(
        DALEX_COMPARE_DIR,
        expected_outputs
      )
    )
)


print(
  output_audit
)


if (
  !all(
    output_audit$Exists
  )
) {

  warning(
    "One or more expected outputs were not created."
  )
}


data.table::fwrite(
  output_audit,
  file.path(
    DALEX_COMPARE_DIR,
    "RA_DALEX_ALE_output_audit.csv"
  )
)


# ======================================================================
# 31. FINAL SUMMARY
# ======================================================================

cat(
  "\n\n============================================================\n",
  "RA LIGHTGBM vs TABPFN DALEX + ALE COMPLETE\n",
  "============================================================\n",
  "Population: original diagnostic FINAL ICE\n",
  "Samples: ",
  nrow(
    dalex_compare_x
  ),
  "\nGenes available to both models: ",
  ncol(
    dalex_compare_x
  ),
  "\n",
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
  "Importance Spearman rho: ",
  sprintf(
    "%.3f",
    importance_spearman
  ),
  "\n",
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
  "\n",
  "============================================================\n",
  sep = ""
)
