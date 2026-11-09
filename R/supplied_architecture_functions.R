# TRACE-Omics | Configurable matched-architecture helpers
# Refactored from the excerpts provided for R/supplied_architecture_functions.R.
# This file may NOT represent every function in the private repository original;
# review the diff against the full original before replacing it.
# R/Python numerical equivalence requires version-pinned environment validation.

foundation_class_weights <- function(outcome, negative_class = "control",
                                     positive_class = "case") {
  y <- factor(as.character(outcome), levels = c(negative_class, positive_class))
  if (anyNA(y)) stop("Unknown outcome label in foundation model.")
  counts <- table(y)
  if (any(counts == 0L)) stop("Both outcome classes are required for class weighting.")
  weights <- length(y) / (2 * as.numeric(counts))
  names(weights) <- names(counts)
  weights
}

select_foundation_rf_features <- function(train_matrix, train_metadata, top_n,
                                          seed = 7832L, rf_prefilter = 3000L,
                                          importance_trees = 250L, threads = 0L,
                                          positive_class = "case", negative_class = "control") {
  trace_require(c("matrixStats", "ranger"))
  if (!is.matrix(train_matrix) || is.null(rownames(train_matrix)) ||
      ncol(train_matrix) != nrow(train_metadata))
    stop("Training matrix must be genes-by-samples and aligned with metadata.")
  if (!"outcome" %in% names(train_metadata)) stop("Metadata require outcome.")
  top_n <- as.integer(top_n)
  rf_prefilter <- as.integer(rf_prefilter)
  if (top_n < 1L || rf_prefilter < top_n || importance_trees < 1L ||
      length(threads) != 1L || !is.finite(threads) || threads < 0L)
    stop("Invalid RF feature selection settings.")
  gene_variance <- matrixStats::rowVars(train_matrix)
  valid <- is.finite(gene_variance) & gene_variance > 0
  train_matrix <- train_matrix[valid, , drop = FALSE]
  gene_variance <- gene_variance[valid]
  if (nrow(train_matrix) > rf_prefilter) {
    ranked <- order(-gene_variance, rownames(train_matrix))
    train_matrix <- train_matrix[head(ranked, rf_prefilter), , drop = FALSE]
  }
  if (nrow(train_matrix) < top_n) stop("Fewer variable genes than requested top_n.")
  x <- as.data.frame(t(train_matrix), check.names = FALSE)
  original_genes <- colnames(x)
  names(x) <- make.names(original_genes, unique = TRUE)
  y <- factor(as.character(train_metadata$outcome),
              levels = c(negative_class, positive_class))
  weights <- foundation_class_weights(y, negative_class, positive_class)
  trace_rng(seed)
  model <- ranger::ranger(x = x, y = y, probability = TRUE,
      importance = "permutation", num.trees = as.integer(importance_trees),
      mtry = max(1L, floor(sqrt(ncol(x)))), class.weights = weights,
      num.threads = as.integer(threads), seed = as.integer(seed))
  values <- model$variable.importance
  values[!is.finite(values)] <- -Inf
  ranking <- order(-values, names(values))
  selected <- match(names(values)[head(ranking, top_n)], names(x))
  if (length(selected) != top_n || anyNA(selected))
    stop("RF feature selection did not return the requested number of genes.")
  original_genes[selected]
}

# Lazy imports prevent source() or parse() from depending on Python packages.
trace_python_module <- function(name, module = NULL) {
  if (!is.null(module)) return(module)
  trace_require("reticulate")
  reticulate::import(name, delay_load = FALSE)
}

to_numpy_x <- function(x, np = NULL) {
  np <- trace_python_module("numpy", np)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)))
    stop("Model predictors must be a finite numeric matrix.")
  np$array(x, dtype = "float32")
}

to_numpy_y <- function(y, np = NULL, positive_class = "case") {
  np <- trace_python_module("numpy", np)
  y <- as.character(y)
  if (!all(y %in% c("control", positive_class)))
    stop("Model outcomes must be control and positive class.")
  np$array(as.integer(y == positive_class), dtype = "int64")
}

fit_lightgbm_object <- function(x_train, y_train, parameters, seed = 7832L,
                                lightgbm_py = NULL, np = NULL, threads = -1L) {
  module <- trace_python_module("lightgbm", lightgbm_py)
  np <- trace_python_module("numpy", np)
  required <- c("n_estimators", "learning_rate", "num_leaves", "max_depth",
                "min_child_samples", "subsample", "colsample_bytree", "reg_alpha", "reg_lambda")
  absent <- setdiff(required, names(parameters))
  if (length(absent)) stop("LightGBM parameter(s) missing: ", paste(absent, collapse = ", "))
  item <- function(k) parameters[[k]][[1L]]
  model <- module$LGBMClassifier(
    objective = "binary",
    n_estimators = as.integer(item("n_estimators")),
    learning_rate = as.numeric(item("learning_rate")),
    num_leaves = as.integer(item("num_leaves")),
    max_depth = as.integer(item("max_depth")),
    min_child_samples = as.integer(item("min_child_samples")),
    subsample = as.numeric(item("subsample")),
    subsample_freq = 1L,
    colsample_bytree = as.numeric(item("colsample_bytree")),
    reg_alpha = as.numeric(item("reg_alpha")),
    reg_lambda = as.numeric(item("reg_lambda")),
    class_weight = "balanced",
    random_state = as.integer(seed),
    n_jobs = as.integer(threads),
    verbosity = -1L
  )
  model$fit(to_numpy_x(x_train, np), to_numpy_y(y_train, np))
  model
}

predict_python_classifier <- function(model, x_test, np = NULL) {
  trace_require("reticulate")
  np <- trace_python_module("numpy", np)
  probabilities <- reticulate::py_to_r(model$predict_proba(to_numpy_x(x_test, np)))
  if (!is.matrix(probabilities)) probabilities <- as.matrix(probabilities)
  classes <- reticulate::py_to_r(model$classes_)
  positive_col <- match(1L, as.integer(classes))
  if (ncol(probabilities) != length(classes) || is.na(positive_col))
    stop("Could not identify probability column for positive class (1).")
  positive <- as.numeric(probabilities[, positive_col])
  if (anyNA(positive) || any(!is.finite(positive)) || any(positive < 0 | positive > 1))
    stop("Python model returned invalid class probabilities.")
  positive
}

# Dissertation TabPFN models were fitted on CPU; device stays configurable.
fit_tabpfn_object <- function(x_train, y_train, seed = 7832L,
                              device = "cpu", tabpfn_py = NULL, np = NULL) {
  module <- trace_python_module("tabpfn", tabpfn_py)
  np <- trace_python_module("numpy", np)
  seed_integer <- as.integer(as.double(seed) %% .Machine$integer.max)
  np$random$seed(seed_integer)
  model <- module$TabPFNClassifier(device = device, seed = seed_integer)
  model$fit(to_numpy_x(x_train, np), to_numpy_y(y_train, np))
  model
}

# Reference architecture grid from the supplied dissertation excerpt.
# Returns a data.frame to avoid loading data.table when merely sourcing code.
trace_lightgbm_grid <- function() {
  data.frame(
    Config_ID = sprintf("LGBM_%02d", seq_len(10L)),
    n_estimators = c(100L, 200L, 200L, 300L, 300L, 150L, 250L, 400L, 250L, 500L),
    learning_rate = c(0.05, 0.03, 0.05, 0.03, 0.02, 0.10, 0.05, 0.02, 0.03, 0.015),
    num_leaves = c(7L, 7L, 15L, 15L, 15L, 7L, 31L, 15L, 31L, 15L),
    max_depth = c(3L, 3L, 4L, 4L, 5L, 3L, 5L, 5L, 6L, 4L),
    min_child_samples = c(10L, 10L, 10L, 15L, 20L, 20L, 15L, 20L, 20L, 10L),
    subsample = c(0.80, 0.80, 0.80, 0.90, 0.80, 0.80, 0.80, 1.00, 0.80, 0.90),
    colsample_bytree = c(0.80, 0.80, 0.80, 0.90, 0.80, 0.80, 0.80, 1.00, 0.90, 0.80),
    reg_alpha = c(0, 0, 0, 0, 0.5, 0, 0.5, 0, 1, 0.5),
    reg_lambda = c(1, 1, 1, 2, 2, 1, 2, 5, 5, 2),
    stringsAsFactors = FALSE
  )
}

lightgbm_grid <- data.table::as.data.table(
  trace_lightgbm_grid()
)

