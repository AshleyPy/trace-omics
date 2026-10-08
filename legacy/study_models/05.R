# ======================================================================
# PART 1
# TB FINAL PUBLICATION PIPELINE
#
# Run AFTER:
#   - 25/50/100 gene optimisation
#   - class-weighted/unweighted optimisation
#   - winner has been locked
#   - FINAL_RF has been fitted
#
# PRIMARY VALIDATION:
#   Strict leave-one-cohort-out predictions
#
# ROC DISPLAY:
#   EMPIRICAL STEPWISE ROC
#
# ROC AUROC CI:
#   DeLong 95% CI
#
# ROC CURVE CI:
#   Pointwise stratified-bootstrap 95% CI
#
# INTERPRETATION:
#   DALEX + repeated Delta-AUROC permutation importance
# ======================================================================


# ======================================================================
# 1. PACKAGES
# ======================================================================

publication_packages <- c(
  "data.table",
  "ggplot2",
  "pROC",
  "PRROC",
  "RColorBrewer",
  "patchwork",
  "scales",
  "knitr",
  "DALEX",
  "ingredients"
)

missing_publication_packages <- publication_packages[
  !vapply(
    publication_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_publication_packages) > 0L) {

  install.packages(
    missing_publication_packages
  )
}

suppressPackageStartupMessages({

  library(data.table)
  library(ggplot2)
  library(pROC)
  library(PRROC)
  library(RColorBrewer)
  library(patchwork)
  library(scales)
  library(knitr)
  library(DALEX)
  library(ingredients)
})


# ======================================================================
# 2. REQUIRED OBJECT AUDIT
# ======================================================================

required_objects <- c(
  "LOCKED_CONFIG_ID",
  "LOCKED_TOP_N",
  "LOCKED_WEIGHTING",
  "LOCKED_THRESHOLD",
  "locked_configuration",
  "locked_loco_predictions",
  "optimisation_metrics",
  "optimisation_summary",
  "diagnostic_ice_predictions",
  "FINAL_RF",
  "FINAL_GENES",
  "final_x_ice",
  "final_ice_counts",
  "final_ice_metadata",
  "development_counts",
  "development_metadata",
  "FINAL_DIR",
  "TABLE_DIR",
  "FIGURE_DIR",
  "MODEL_DIR",
  "SET1",
  "POSITIVE_CLASS",
  "NEGATIVE_CLASS"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    FUN.VALUE = logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0L) {

  stop(
    "Missing required objects:\n",
    paste(
      missing_objects,
      collapse = "\n"
    )
  )
}


# ======================================================================
# 3. OUTPUT DIRECTORIES
# ======================================================================

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

DALEX_DIR <- file.path(
  FINAL_DIR,
  "DALEX"
)

dir.create(
  DALEX_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 4. SETTINGS
# ======================================================================

SEED_PUBLICATION <- 7832L

ROC_BOOTSTRAPS <- 2000L

ROC_SPECIFICITY_GRID <- seq(
  0,
  1,
  length.out = 101L
)

B_IMPORTANCE <- 100L

TOP_IMPORTANCE_GENES <- min(
  15L,
  length(
    FINAL_GENES
  )
)

PROFILE_N <- min(
  6L,
  length(
    FINAL_GENES
  )
)


# ======================================================================
# 5. COHORT ORDER
# ======================================================================

cohort_order <- c(
  "ACS",
  "Catalysis",
  "GC6"
)

cohort_order <- cohort_order[
  cohort_order %in%
    unique(
      locked_loco_predictions$Target_cohort
    )
]

if (length(cohort_order) != 3L) {

  stop(
    "Expected ACS, Catalysis and GC6 in strict LOCO predictions."
  )
}


# ======================================================================
# 6. SET1 COLOURS
# ======================================================================

cohort_colours <- stats::setNames(
  SET1[
    seq_along(
      cohort_order
    )
  ],
  cohort_order
)
# ======================================================================
# PART 2
# PUBLICATION THEME
# ======================================================================

publication_theme <- function() {

  ggplot2::theme_bw(
    base_size = 11
  ) +

    ggplot2::theme(

      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13
      ),

      plot.subtitle = ggplot2::element_text(
        size = 9.5
      ),

      axis.title = ggplot2::element_text(
        face = "bold"
      ),

      legend.title = ggplot2::element_text(
        face = "bold"
      ),

      legend.position = "right",

      panel.grid.minor = ggplot2::element_blank(),

      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black"
      ),

      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )
}
# ======================================================================
# PART 3A
# THRESHOLD METRICS
# ======================================================================

calculate_threshold_metrics_publication <- function(
  truth,
  probability,
  threshold
) {

  truth <- as.character(
    truth
  )

  probability <- as.numeric(
    probability
  )

  response <- ifelse(
    probability >= threshold,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  )

  sensitivity <- mean(
    response[
      truth == POSITIVE_CLASS
    ] ==
      POSITIVE_CLASS
  )

  specificity <- mean(
    response[
      truth == NEGATIVE_CLASS
    ] ==
      NEGATIVE_CLASS
  )

  data.table::data.table(

    BACC = mean(
      c(
        sensitivity,
        specificity
      )
    ),

    Sensitivity = sensitivity,

    Specificity = specificity,

    Brier = mean(
      (
        probability -
          as.integer(
            truth == POSITIVE_CLASS
          )
      )^2
    )
  )
}


# ======================================================================
# PART 3B
# CALIBRATION BINS
# ======================================================================

make_calibration_data <- function(
  truth,
  probability,
  label,
  number_of_bins = 6L
) {

  calibration_dt <- data.table::data.table(

    truth = as.integer(
      truth == POSITIVE_CLASS
    ),

    probability = as.numeric(
      probability
    )
  )

  calibration_dt[
    ,
    Probability_rank := data.table::frank(
      probability,
      ties.method = "average"
    )
  ]

  calibration_dt[
    ,
    Bin := ceiling(
      Probability_rank /
        .N *
        number_of_bins
    )
  ]

  calibration_dt[
    Bin < 1L,
    Bin := 1L
  ]

  calibration_dt[
    Bin > number_of_bins,
    Bin := number_of_bins
  ]

  result <- calibration_dt[
    ,
    .(
      n = .N,
      Cases = sum(
        truth
      ),
      Mean_predicted = mean(
        probability
      ),
      Observed = mean(
        truth
      )
    ),
    by = Bin
  ]

  result[
    ,
    Lower95 := vapply(
      seq_len(.N),
      function(row_index) {

        stats::binom.test(
          x = Cases[[row_index]],
          n = n[[row_index]]
        )$conf.int[[1L]]
      },
      FUN.VALUE = numeric(1)
    )
  ]

  result[
    ,
    Upper95 := vapply(
      seq_len(.N),
      function(row_index) {

        stats::binom.test(
          x = Cases[[row_index]],
          n = n[[row_index]]
        )$conf.int[[2L]]
      },
      FUN.VALUE = numeric(1)
    )
  ]

  result[
    ,
    Group := label
  ]

  result
}


# ======================================================================
# PART 3C
# CALIBRATION INTERCEPT / SLOPE
# ======================================================================

calculate_calibration_statistics <- function(
  truth,
  probability
) {

  truth_numeric <- as.integer(
    truth == POSITIVE_CLASS
  )

  probability_clamped <- pmin(
    pmax(
      as.numeric(
        probability
      ),
      1e-6
    ),
    1 - 1e-6
  )

  linear_predictor <- stats::qlogis(
    probability_clamped
  )

  slope_model <- stats::glm(
    truth_numeric ~ linear_predictor,
    family = stats::binomial()
  )

  intercept_model <- stats::glm(
    truth_numeric ~ offset(
      linear_predictor
    ),
    family = stats::binomial()
  )

  data.table::data.table(

    Calibration_intercept = as.numeric(
      stats::coef(
        intercept_model
      )[[1L]]
    ),

    Calibration_slope = as.numeric(
      stats::coef(
        slope_model
      )[["linear_predictor"]]
    ),

    Brier = mean(
      (
        truth_numeric -
          probability_clamped
      )^2
    )
  )
}
