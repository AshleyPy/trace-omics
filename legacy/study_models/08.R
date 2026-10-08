# ======================================================================
# PR BOOTSTRAP HELPER
#
# Pointwise 95% CI for precision at fixed recall.
#
# Bootstrap is stratified:
#   cases are resampled from cases
#   controls are resampled from controls
#
# This preserves cohort case prevalence in every bootstrap replicate.
# ======================================================================

PR_BOOTSTRAPS <- 2000L

PR_RECALL_GRID <- seq(
  0,
  1,
  length.out = 101L
)


get_pr_precision_on_grid <- function(
  truth,
  probability,
  recall_grid
) {

  truth <- as.character(
    truth
  )

  current_pr <- PRROC::pr.curve(

    scores.class0 =
      probability[
        truth ==
          POSITIVE_CLASS
      ],

    scores.class1 =
      probability[
        truth ==
          NEGATIVE_CLASS
      ],

    curve =
      TRUE
  )


  current_curve <- data.table::data.table(

    Recall =
      as.numeric(
        current_pr$curve[
          ,
          1L
        ]
      ),

    Precision =
      as.numeric(
        current_pr$curve[
          ,
          2L
        ]
      )
  )


  current_curve <- current_curve[
    is.finite(
      Recall
    ) &
      is.finite(
        Precision
      )
  ]


  # Multiple PR points can have the same recall.
  # Retain the largest precision at each recall value.
  current_curve <- current_curve[
    ,
    .(
      Precision =
        max(
          Precision,
          na.rm = TRUE
        )
    ),
    by =
      Recall
  ]


  data.table::setorder(
    current_curve,
    Recall
  )


  if (
    nrow(
      current_curve
    ) <
      2L
  ) {

    return(
      rep(
        NA_real_,
        length(
          recall_grid
        )
      )
    )
  }


  stats::approx(

    x =
      current_curve$Recall,

    y =
      current_curve$Precision,

    xout =
      recall_grid,

    method =
      "linear",

    rule =
      2,

    ties =
      "ordered"
  )$y
}


bootstrap_pr_ci <- function(
  truth,
  probability,
  recall_grid = PR_RECALL_GRID,
  bootstraps = PR_BOOTSTRAPS,
  seed = 7832L
) {

  truth <- as.character(
    truth
  )

  probability <- as.numeric(
    probability
  )


  case_index <- which(
    truth ==
      POSITIVE_CLASS
  )


  control_index <- which(
    truth ==
      NEGATIVE_CLASS
  )


  if (
    length(
      case_index
    ) <
      2L ||
      length(
        control_index
      ) <
        2L
  ) {

    stop(
      "PR bootstrap requires at least two cases and two controls."
    )
  }


  bootstrap_matrix <- matrix(
    NA_real_,
    nrow = bootstraps,
    ncol = length(
      recall_grid
    )
  )


  set.seed(
    seed
  )


  for (
    bootstrap_index in seq_len(
      bootstraps
    )
  ) {

    sampled_case_index <- sample(
      case_index,
      size = length(
        case_index
      ),
      replace = TRUE
    )


    sampled_control_index <- sample(
      control_index,
      size = length(
        control_index
      ),
      replace = TRUE
    )


    sampled_index <- c(
      sampled_case_index,
      sampled_control_index
    )


    bootstrap_matrix[
      bootstrap_index,
    ] <- get_pr_precision_on_grid(

      truth =
        truth[
          sampled_index
        ],

      probability =
        probability[
          sampled_index
        ],

      recall_grid =
        recall_grid
    )
  }


  data.table::data.table(

    Recall =
      recall_grid,

    Precision_lower =
      apply(
        bootstrap_matrix,
        2L,
        stats::quantile,
        probs = 0.025,
        na.rm = TRUE
      ),

    Precision_median =
      apply(
        bootstrap_matrix,
        2L,
        stats::quantile,
        probs = 0.500,
        na.rm = TRUE
      ),

    Precision_upper =
      apply(
        bootstrap_matrix,
        2L,
        stats::quantile,
        probs = 0.975,
        na.rm = TRUE
      )
  )
}
