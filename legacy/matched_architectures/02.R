# ======================================================================
# RA LIGHTGBM vs TABPFN
# DALEX PERMUTATION IMPORTANCE + ALE
#
# USES ONLY THE ALREADY-COMPLETED DIAGNOSTIC FINAL ICE TEST
#
# NO:
#   - preprocessing
#   - feature selection
#   - tuning
#   - refitting
#   - nested CV
#   - LOCO
#   - diagnostic ICE re-testing
# ======================================================================


# ======================================================================
# 0. PATHS
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
# 1. LOCK CORRECT RA PYTHON ENVIRONMENT
#
# RUN IN A FRESH R SESSION.
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
    "Could not find Python for environment:\n",
    FOUNDATION_CONDA_ENV
  )
}


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


  requested_python <- normalizePath(
    FOUNDATION_PYTHON,
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
      "Wrong Python is already active.\n\n",
      "Restart RStudio with Ctrl + Shift + F10, then run this script ",
      "from the beginning.\n\n",
      "Active:\n",
      active_python,
      "\n\nRequired:\n",
      requested_python
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


# ======================================================================
# 2. PACKAGES
# ======================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "RColorBrewer",
  "pROC",
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
  ) >
    0L
) {

  stop(
    "Install these packages first:\n",
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
  library(pROC)
  library(reticulate)
  library(DALEX)
  library(ingredients)
  library(patchwork)

})


np <- reticulate::import(
  "numpy",
  convert = FALSE
)


# ======================================================================
# 3. LOAD SAVED RA EXPERIMENT STATE
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


cat(
  "\n============================================================\n",
  "RA FOUNDATION STATE LOADED\n",
  "============================================================\n",
  "Disease: ",
  foundation_state$disease,
  "\nRepresentation: ",
  foundation_state$representation,
  "\nSelected genes: ",
  length(
    ICE_SELECTED_GENES
  ),
  "\nPositive class: ",
  POSITIVE_CLASS,
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 4. LOAD EXACT OLD DIAGNOSTIC ICE MATRIX
#
# This is the matrix that was already passed to LightGBM and TabPFN
# during the completed diagnostic FINAL ICE test.
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
    "Could not find the saved diagnostic ICE matrix:\n",
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
    "The saved diagnostic ICE matrix genes do not match ",
    "the saved RA experiment state."
  )
}


cat(
  "\n============================================================\n",
  "OLD RA DIAGNOSTIC ICE MATRIX LOADED\n",
  "============================================================\n",
  "Samples: ",
  nrow(
    ice_x_test
  ),
  "\nGenes: ",
  ncol(
    ice_x_test
  ),
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 5. LOAD OLD DIAGNOSTIC ICE PREDICTIONS
#
# These provide the original sample IDs, truth labels and probabilities.
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


required_prediction_columns <- c(
  "Model",
  "sample_id",
  "truth",
  "probability_case"
)


missing_prediction_columns <- setdiff(
  required_prediction_columns,
  names(
    ice_predictions
  )
)


if (
  length(
    missing_prediction_columns
  ) >
    0L
) {

  stop(
    "Prediction file is missing:\n",
    paste(
      missing_prediction_columns,
      collapse = ", "
    )
  )
}


# ======================================================================
# 6. ALIGN OLD TRUTH LABELS TO THE SAVED ICE MATRIX
# ======================================================================

truth_dt <- unique(
  ice_predictions[
    ,
    .(
      sample_id,
      truth
    )
  ]
)


if (
  anyDuplicated(
    truth_dt$sample_id
  ) >
    0L
) {

  stop(
    "A sample has inconsistent truth labels."
  )
}


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
    "Could not align the saved ICE matrix to the old ",
    "diagnostic-ICE predictions."
  )
}


ice_truth <- truth_dt$truth[
  sample_match
]


dalex_compare_y <- as.integer(
  ice_truth ==
    POSITIVE_CLASS
)


dalex_compare_x <- as.data.frame(
  ice_x_test,
  check.names = FALSE
)


# ======================================================================
# 7. LOAD THE ALREADY-FITTED FINAL MODELS
#
# NO FITTING OCCURS HERE.
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


message(
  "PASS: saved RA LightGBM and TabPFN models loaded."
)


# ======================================================================
# 8. PREDICTION FUNCTIONS
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


  probability_py <- model$predict_proba(
    newdata_py
  )


  probability_matrix <- reticulate::py_to_r(
    probability_py
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


  probability_py <- model$predict_proba(
    newdata_py
  )


  probability_matrix <- reticulate::py_to_r(
    probability_py
  )


  as.numeric(
    probability_matrix[, 2L]
  )
}


# ======================================================================
# 9. AUDIT AGAINST OLD DIAGNOSTIC ICE PROBABILITIES
#
# This does NOT create a new test.
# It only confirms that the reloaded saved model reproduces the old result.
# ======================================================================

audit_model_predictions <- function(
  model,
  model_name,
  predict_function
) {

  original_dt <- ice_predictions[
    Model ==
      model_name
  ]


  original_match <- match(
    rownames(
      ice_x_test
    ),
    original_dt$sample_id
  )


  if (
    anyNA(
      original_match
    )
  ) {

    stop(
      "Could not align old ",
      model_name,
      " predictions."
    )
  }


  original_probability <- original_dt$probability_case[
    original_match
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
          restored_probability -
            original_probability
        )
      ),

    Mean_absolute_difference =
      mean(
        abs(
          restored_probability -
            original_probability
        )
      )
  )
}


prediction_audit <- data.table::rbindlist(
  list(

    audit_model_predictions(
      FINAL_LGBM_MODEL,
      "LightGBM",
      predict_lightgbm_dalex
    ),

    audit_model_predictions(
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
    "Reloaded model probabilities do not reproduce the old ",
    "diagnostic ICE results. DALEX has been stopped."
  )
}


message(
  "PASS: restored models reproduce the old RA diagnostic ICE predictions."
)


data.table::fwrite(
  prediction_audit,
  file.path(
    DALEX_COMPARE_DIR,
    "00_RA_saved_model_prediction_audit.csv"
  )
)


# ======================================================================
# 10. BUILD DALEX EXPLAINERS ON OLD DIAGNOSTIC ICE
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
# 11. SETTINGS
# ======================================================================

SEED_FOUNDATION <- 7832L

IMPORTANCE_PERMUTATIONS <- 100L

TOP_COMPARISON_GENES <- 20L

PROFILE_N <- 6L


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


# ======================================================================
# 12. PUBLICATION THEME
# ======================================================================

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
          face = "bold"
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
          face = "bold"
        )
    )
}


# ======================================================================
# 13. REPEATED DELTA-AUROC PERMUTATION IMPORTANCE
#
# Same method as the previous LightGBM vs TabPFN comparison:
#
# Delta AUROC = original AUROC - AUROC after permuting one gene
#
# 100 repeated permutations per gene.
# ======================================================================

calculate_model_permutation_importance <- function(
  model,
  model_name,
  predict_function,
  x,
  truth,
  permutations,
  seed
) {

  x <- as.data.frame(
    x,
    check.names = FALSE
  )


  truth_factor <- factor(
    truth,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )


  baseline_probability <- predict_function(
    model,
    x
  )


  baseline_roc <- pROC::roc(
    response = truth_factor,
    predictor = baseline_probability,
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    ),
    direction = "<",
    quiet = TRUE
  )


  baseline_auc <- as.numeric(
    pROC::auc(
      baseline_roc
    )
  )


  raw_rows <- vector(
    mode = "list",
    length = ncol(
      x
    )
  )


  for (
    gene_index in seq_len(
      ncol(
        x
      )
    )
  ) {

    current_gene <- colnames(
      x
    )[[gene_index]]


    original_values <- x[
      ,
      gene_index
    ]


    gene_rows <- vector(
      mode = "list",
      length = permutations
    )


    for (
      permutation_index in seq_len(
        permutations
      )
    ) {

      set.seed(
        seed +
          gene_index *
            10000L +
          permutation_index
      )


      x_permuted <- x


      x_permuted[
        ,
        gene_index
      ] <- sample(
        original_values,
        size = length(
          original_values
        ),
        replace = FALSE
      )


      permuted_probability <- predict_function(
        model,
        x_permuted
      )


      permuted_roc <- pROC::roc(
        response = truth_factor,
        predictor = permuted_probability,
        levels = c(
          NEGATIVE_CLASS,
          POSITIVE_CLASS
        ),
        direction = "<",
        quiet = TRUE
      )


      permuted_auc <- as.numeric(
        pROC::auc(
          permuted_roc
        )
      )


      gene_rows[[permutation_index]] <-
        data.table::data.table(

          Model =
            model_name,

          Gene =
            current_gene,

          Permutation =
            permutation_index,

          Baseline_AUROC =
            baseline_auc,

          Permuted_AUROC =
            permuted_auc,

          Delta_AUROC =
            baseline_auc -
            permuted_auc
        )
    }


    raw_rows[[gene_index]] <- data.table::rbindlist(
      gene_rows
    )


    if (
      gene_index %% 5L ==
        0L ||
        gene_index ==
          ncol(
            x
          )
    ) {

      message(
        model_name,
        ": ",
        gene_index,
        "/",
        ncol(
          x
        ),
        " genes complete."
      )
    }
  }


  raw_dt <- data.table::rbindlist(
    raw_rows
  )


  summary_dt <- raw_dt[
    ,
    .(

      Baseline_AUROC =
        unique(
          Baseline_AUROC
        )[[1L]],

      Mean_permuted_AUROC =
        mean(
          Permuted_AUROC,
          na.rm = TRUE
        ),

      Mean_Delta_AUROC =
        mean(
          Delta_AUROC,
          na.rm = TRUE
        ),

      SD_Delta_AUROC =
        stats::sd(
          Delta_AUROC,
          na.rm = TRUE
        ),

      Lower95_Delta_AUROC =
        as.numeric(
          stats::quantile(
            Delta_AUROC,
            probs = 0.025,
            na.rm = TRUE
          )
        ),

      Upper95_Delta_AUROC =
        as.numeric(
          stats::quantile(
            Delta_AUROC,
            probs = 0.975,
            na.rm = TRUE
          )
        ),

      Positive_permutation_fraction =
        mean(
          Delta_AUROC >
            0,
          na.rm = TRUE
        )
    ),
    by = .(
      Model,
      Gene
    )
  ]


  data.table::setorder(
    summary_dt,
    Model,
    -Mean_Delta_AUROC,
    Gene
  )


  summary_dt[
    ,
    Rank :=
      seq_len(
        .N
      ),
    by =
      Model
  ]


  list(
    raw = raw_dt,
    summary = summary_dt
  )
}


# ======================================================================
# 14. LIGHTGBM DALEX IMPORTANCE
# ======================================================================

cat(
  "\n============================================================\n",
  "RA LIGHTGBM DALEX IMPORTANCE\n",
  "============================================================\n"
)


lgbm_importance_result <- calculate_model_permutation_importance(

  model =
    FINAL_LGBM_MODEL,

  model_name =
    "LightGBM",

  predict_function =
    predict_lightgbm_dalex,

  x =
    dalex_compare_x,

  truth =
    ice_truth,

  permutations =
    IMPORTANCE_PERMUTATIONS,

  seed =
    SEED_FOUNDATION +
      88200001L
)


lgbm_importance_summary <-
  lgbm_importance_result$summary


# ======================================================================
# 15. TABPFN DALEX IMPORTANCE
# ======================================================================

cat(
  "\n============================================================\n",
  "RA TABPFN DALEX IMPORTANCE\n",
  "============================================================\n"
)


tabpfn_importance_result <- calculate_model_permutation_importance(

  model =
    FINAL_TABPFN_MODEL,

  model_name =
    "TabPFN",

  predict_function =
    predict_tabpfn_dalex,

  x =
    dalex_compare_x,

  truth =
    ice_truth,

  permutations =
    IMPORTANCE_PERMUTATIONS,

  seed =
    SEED_FOUNDATION +
      88200002L
)


tabpfn_importance_summary <-
  tabpfn_importance_result$summary


# ======================================================================
# 16. SAVE IMPORTANCE RESULTS
# ======================================================================

comparison_importance_raw <- data.table::rbindlist(
  list(
    lgbm_importance_result$raw,
    tabpfn_importance_result$raw
  ),
  use.names = TRUE,
  fill = TRUE
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
  comparison_importance_raw,
  file.path(
    DALEX_COMPARE_DIR,
    "01_RA_LightGBM_TabPFN_permutation_importance_raw.csv"
  )
)


data.table::fwrite(
  comparison_importance_summary,
  file.path(
    DALEX_COMPARE_DIR,
    "02_RA_LightGBM_TabPFN_permutation_importance_summary.csv"
  )
)


# ======================================================================
# 17. IMPORTANCE PLOT FUNCTION
# ======================================================================

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
      title = paste0(
        model_name,
        " global feature importance"
      ),
      subtitle = paste0(
        "Top ",
        top_n,
        " genes ranked by repeated permutation \u0394AUROC"
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
      legend.position = "top",
      legend.justification = "left",
      axis.text.y = ggplot2::element_text(
        size = 9,
        face = "italic",
        colour = "black"
      ),
      panel.grid.major.y =
        ggplot2::element_blank()
    )
}


# ======================================================================
# 18. LIGHTGBM + TABPFN IMPORTANCE FIGURES
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
      "Predictive feature reliance in RA: LightGBM versus TabPFN",

    subtitle =
      paste0(
        "Post-hoc importance in the identical diagnostic held-out population; ",
        IMPORTANCE_PERMUTATIONS,
        " permutations per gene"
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


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.pdf"
  ),
  plot = p_model_importance_comparison,
  width = 13.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  device = grDevices::cairo_pdf
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_DALEX_LightGBM_vs_TabPFN_importance.png"
  ),
  plot = p_model_importance_comparison,
  width = 13.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 19. TOP SIX GENES FOR ALE
#
# Same approach as the previous comparison:
# each model gets its own six highest-ranked genes.
# ======================================================================

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
  "\nLightGBM ALE genes:\n"
)


print(
  LGBM_PROFILE_GENES
)


cat(
  "\nTabPFN ALE genes:\n"
)


print(
  TABPFN_PROFILE_GENES
)


# ======================================================================
# 20. LIGHTGBM ALE
# ======================================================================

set.seed(
  SEED_FOUNDATION
)


lgbm_ale <- ingredients::accumulated_dependency(

  lgbm_explainer,

  variables =
    LGBM_PROFILE_GENES,

  N =
    NULL
)


# ======================================================================
# 21. TABPFN ALE
# ======================================================================

set.seed(
  SEED_FOUNDATION
)


tabpfn_ale <- ingredients::accumulated_dependency(

  tabpfn_explainer,

  variables =
    TABPFN_PROFILE_GENES,

  N =
    NULL
)


# ======================================================================
# 22. SAVE ALE RESULTS
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


data.table::fwrite(
  data.table::as.data.table(
    lgbm_ale
  ),
  file.path(
    DALEX_COMPARE_DIR,
    "RA_LightGBM_ALE.csv"
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
# 23. ALE FIGURES
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
      "Gene expression",
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


p_tabpfn_ale <- plot(
  tabpfn_ale
) +

  ggplot2::labs(
    title =
      "TabPFN expression-response profiles",
    subtitle =
      "Accumulated local effects for the six most influential RA genes",
    x =
      "Gene expression",
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


# ======================================================================
# 24. COMBINED ALE FIGURE
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
      "RA expression-response profiles: LightGBM versus TabPFN",
    subtitle =
      "Accumulated local effects in the identical diagnostic held-out population",
    tag_levels =
      "A"
  )


print(
  p_ale_comparison
)


ggplot2::ggsave(
  filename = file.path(
    DALEX_COMPARE_DIR,
    "Figure_RA_ALE_LightGBM_vs_TabPFN.pdf"
  ),
  plot = p_ale_comparison,
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
  plot = p_ale_comparison,
  width = 13.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ======================================================================
# 25. FINAL AUDIT
# ======================================================================

cat(
  "\n\n============================================================\n",
  "RA DALEX + ALE COMPLETE\n",
  "============================================================\n",
  "Interpretation population: ORIGINAL diagnostic FINAL ICE\n",
  "Samples: ",
  nrow(
    dalex_compare_x
  ),
  "\nGenes: ",
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
  "Permutation repeats per gene: ",
  IMPORTANCE_PERMUTATIONS,
  "\n\n",
  "NO preprocessing rerun\n",
  "NO feature selection rerun\n",
  "NO model tuning rerun\n",
  "NO model fitting rerun\n",
  "NO diagnostic ICE test rerun\n\n",
  "Saved to:\n",
  DALEX_COMPARE_DIR,
  "\n============================================================\n",
  sep = ""
)

