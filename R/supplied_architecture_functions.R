# Supplied matched architecture functions and ten-row grid.
foundation_class_weights <- function(outcome) {

  y <- factor(
    as.character(outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  counts <- table(y)

  weights <- length(y) /
    (2 * as.numeric(counts))

  names(weights) <- names(counts)

  weights
}

select_foundation_rf_features <- function(
  train_matrix,
  train_metadata,
  top_n,
  seed
) {

  gene_variance <- matrixStats::rowVars(
    train_matrix
  )

  valid <- is.finite(gene_variance) &
    gene_variance > 0

  train_matrix <- train_matrix[
    valid,
    ,
    drop = FALSE
  ]

  gene_variance <- gene_variance[
    valid
  ]

  if (nrow(train_matrix) > RF_PREFILTER) {
    variance_order <- order(
      gene_variance,
      decreasing = TRUE
    )

    keep <- variance_order[
      seq_len(RF_PREFILTER)
    ]

    train_matrix <- train_matrix[
      keep,
      ,
      drop = FALSE
    ]
  }

  x <- as.data.frame(
    t(train_matrix),
    check.names = FALSE
  )

  y <- factor(
    as.character(train_metadata$outcome),
    levels = c(
      NEGATIVE_CLASS,
      POSITIVE_CLASS
    )
  )

  class_weights <- foundation_class_weights(y)

  set.seed(seed)

  rf <- ranger::ranger(
    x = x,
    y = y,
    probability = TRUE,
    importance = "permutation",
    num.trees = RF_IMPORTANCE_TREES,
    mtry = max(
      1L,
      floor(sqrt(ncol(x)))
    ),
    class.weights = class_weights,
    num.threads = 0L,
    seed = seed
  )

  importance_values <- rf$variable.importance

  importance_values[
    !is.finite(importance_values)
  ] <- -Inf

  selected_genes <- names(
    head(
      sort(
        importance_values,
        decreasing = TRUE
      ),
      top_n
    )
  )

  if (length(selected_genes) != top_n) {
    stop(
      "RF feature selection returned ",
      length(selected_genes),
      " genes; expected ",
      top_n,
      "."
    )
  }

  selected_genes
}

to_numpy_x <- function(x) {
  np$array(
    x,
    dtype = "float32"
  )
}

to_numpy_y <- function(y) {

  y_numeric <- as.integer(
    as.character(y) ==
      POSITIVE_CLASS
  )

  np$array(
    y_numeric,
    dtype = "int64"
  )
}

fit_lightgbm_object <- function(
  x_train,
  y_train,
  parameters,
  seed
) {

  model <- lightgbm_py$LGBMClassifier(
    objective = "binary",
    n_estimators = as.integer(parameters$n_estimators[[1L]]),
    learning_rate = parameters$learning_rate[[1L]],
    num_leaves = as.integer(parameters$num_leaves[[1L]]),
    max_depth = as.integer(parameters$max_depth[[1L]]),
    min_child_samples = as.integer(parameters$min_child_samples[[1L]]),
    subsample = parameters$subsample[[1L]],
    subsample_freq = 1L,
    colsample_bytree = parameters$colsample_bytree[[1L]],
    reg_alpha = parameters$reg_alpha[[1L]],
    reg_lambda = parameters$reg_lambda[[1L]],
    class_weight = "balanced",
    random_state = as.integer(seed),
    n_jobs = -1L,
    verbosity = -1L
  )

  model$fit(
    to_numpy_x(x_train),
    to_numpy_y(y_train)
  )

  model
}

predict_python_classifier <- function(
  model,
  x_test
) {

  probability_py <- model$predict_proba(
    to_numpy_x(x_test)
  )

  probability_matrix <- reticulate::py_to_r(
    probability_py
  )

  if (!is.matrix(probability_matrix)) {
    probability_matrix <- as.matrix(
      probability_matrix
    )
  }

  if (ncol(probability_matrix) != 2L) {
    stop(
      "Expected two probability columns; received ",
      ncol(probability_matrix),
      "."
    )
  }

  probability_case <- as.numeric(
    probability_matrix[
      ,
      2L
    ]
  )

  if (any(!is.finite(probability_case))) {
    stop("Model returned non-finite probabilities.")
  }

  probability_case
}

fit_tabpfn_object <- function(
  x_train,
  y_train,
  seed
) {

  seed_integer <- as.integer(
    seed %%
      .Machine$integer.max
  )

  np$random$seed(
    seed_integer
  )

  model <- tabpfn_py$TabPFNClassifier(
    device = TABPFN_DEVICE,
    seed = seed_integer
  )

  model$fit(
    to_numpy_x(x_train),
    to_numpy_y(y_train)
  )

  model
}

lightgbm_grid <- data.table::data.table(

  Config_ID = paste0(
    "LGBM_",
    sprintf("%02d", 1:10)
  ),

  n_estimators = c(
    100L,
    200L,
    200L,
    300L,
    300L,
    150L,
    250L,
    400L,
    250L,
    500L
  ),

  learning_rate = c(
    0.05,
    0.03,
    0.05,
    0.03,
    0.02,
    0.10,
    0.05,
    0.02,
    0.03,
    0.015
  ),

  num_leaves = c(
    7L,
    7L,
    15L,
    15L,
    15L,
    7L,
    31L,
    15L,
    31L,
    15L
  ),

  max_depth = c(
    3L,
    3L,
    4L,
    4L,
    5L,
    3L,
    5L,
    5L,
    6L,
    4L
  ),

  min_child_samples = c(
    10L,
    10L,
    10L,
    15L,
    20L,
    20L,
    15L,
    20L,
    20L,
    10L
  ),

  subsample = c(
    0.80,
    0.80,
    0.80,
    0.90,
    0.80,
    0.80,
    0.80,
    1.00,
    0.80,
    0.90
  ),

  colsample_bytree = c(
    0.80,
    0.80,
    0.80,
    0.90,
    0.80,
    0.80,
    0.80,
    1.00,
    0.90,
    0.80
  ),

  reg_alpha = c(
    0,
    0,
    0,
    0,
    0.5,
    0,
    0.5,
    0,
    1,
    0.5
  ),

  reg_lambda = c(
    1,
    1,
    1,
    2,
    2,
    1,
    2,
    5,
    5,
    2
  )
)

print(lightgbm_grid)
