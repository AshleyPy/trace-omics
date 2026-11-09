# TRACE-Omics | Reusable cohort RF and threshold-selection helpers.
# Refactored from the supplied 30-line historical study_models.R.
# The threshold is fitted on DEVELOPMENT out-of-fold predictions only,
# never by inspecting performance on an untouched holdout or LOCO test cohort.

trace_rf <- function(train, new, meta, new_meta, representation, genes,
                     weighted = FALSE, seed = 7832L, parameters = list(),
                     config = list(), project_root = NULL) {
  trace_require(c("ranger", "matrixStats"))
  opt <- function(name, default) if (is.null(parameters[[name]])) default else parameters[[name]]
  prefilter <- as.integer(opt("prefilter", 3000L))
  screening_trees <- as.integer(opt("screening_trees", 250L))
  final_trees <- as.integer(opt("final_trees", 1000L))
  min_node_size <- as.integer(opt("min_node_size", 5L))
  threads <- as.integer(opt("threads", 1L))
  if (genes < 1L || prefilter < genes || screening_trees < 1L || final_trees < 1L ||
      min_node_size < 1L || threads < 1L) stop("Invalid RF configuration.")
  z <- trace_representation(representation, train, new, meta, new_meta,
                            seed = seed, config = config, project_root = project_root)
  v <- matrixStats::rowVars(z$train)
  keep <- which(is.finite(v) & v > 0)
  keep <- keep[order(-v[keep], rownames(z$train)[keep])]
  keep <- head(keep, prefilter)
  if (length(keep) < genes) stop("Fewer variable training genes than requested.")
  x <- as.data.frame(t(z$train[keep, , drop = FALSE]), check.names = FALSE)
  names(x) <- make.names(names(x), unique = TRUE)
  y <- factor(as.character(meta$outcome), levels = c("control", "case"))
  if (anyNA(y) || length(unique(y)) < 2L) stop("Both outcome classes are needed.")
  weights <- if (isTRUE(weighted)) {
    n <- table(y)
    if (any(n == 0L)) stop("Cannot weight an absent class.")
    setNames(length(y) / (2 * as.numeric(n)), names(n))
  } else NULL
  trace_rng(seed)
  screen <- ranger::ranger(x = x, y = y, probability = TRUE,
      num.trees = screening_trees,
      mtry = max(1L, floor(sqrt(ncol(x)))), importance = "permutation",
      class.weights = weights, num.threads = threads, seed = as.integer(seed))
  imp <- screen$variable.importance
  ord <- order(-imp, names(imp), na.last = NA)
  take <- head(ord, as.integer(genes))
  if (length(take) != genes) stop("RF importance produced too few features.")
  selected <- rownames(z$train)[keep[take]]
  xt <- x[, take, drop = FALSE]
  xv <- as.data.frame(t(z$new[selected, , drop = FALSE]), check.names = FALSE)
  names(xv) <- names(xt)
  model <- ranger::ranger(x = xt, y = y, probability = TRUE,
      num.trees = final_trees,
      mtry = max(1L, floor(sqrt(ncol(xt)))), min.node.size = min_node_size,
      class.weights = weights, num.threads = threads, seed = as.integer(seed))
  probability <- as.numeric(predict(model, data = xv)$predictions[, "case"])
  if (length(probability) != nrow(xv) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) stop("Invalid RF case probabilities.")
  list(probability = probability, genes = selected, model = model,
       representation_state = z$state)
}

trace_threshold <- function(predictions, thresholds = seq(0.05, 0.95, 0.01),
                            group_col = "study", truth_col = "outcome",
                            probability_col = "probability") {
  predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
  required <- c(group_col, truth_col, probability_col)
  if (!all(required %in% names(predictions)))
    stop("Development predictions lack: ", paste(setdiff(required, names(predictions)), collapse = ", "))
  if (anyNA(predictions[, required]) ||
      any(predictions[[probability_col]] < 0 | predictions[[probability_col]] > 1))
    stop("Invalid development predictions.")
  if (any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1))
    stop("Thresholds must lie between 0 and 1.")
  grouped <- split(predictions, predictions[[group_col]])
  bacc <- function(y, p, cutoff) {
    y <- as.character(y)
    if (!all(c("control", "case") %in% unique(y)))
      stop("Each development group requires both classes for balanced accuracy.")
    pred <- p >= cutoff
    mean(c(mean(pred[y == "case"]), mean(!pred[y == "control"])))
  }
  scores <- t(vapply(thresholds, function(cutoff) {
    by_group <- vapply(grouped, function(z) bacc(z[[truth_col]], z[[probability_col]], cutoff), numeric(1))
    c(mean = mean(by_group), worst = min(by_group))
  }, numeric(2)))
  thresholds[order(-scores[, "mean"], -scores[, "worst"], abs(thresholds - 0.5))[1L]]
}
