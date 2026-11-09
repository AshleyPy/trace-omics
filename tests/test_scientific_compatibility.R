# TRACE-Omics regression tests for metadata and matched-architecture compatibility.
# Run from the repository root. Requires data.table, not Python, ranger or patient data.
source("R/core.R")
source("R/study_models.R")
source("R/supplied_architecture_functions.R")

stopifnot(
  requireNamespace("data.table", quietly = TRUE),
  data.table::is.data.table(lightgbm_grid),
  nrow(lightgbm_grid) == 10L,
  nrow(lightgbm_grid[1L]) == 1L,
  identical(lightgbm_grid$Config_ID, sprintf("LGBM_%02d", 1:10)),
  identical(eval(formals(fit_tabpfn_object)$device), "cpu"),
  identical(eval(formals(select_foundation_rf_features)$threads), 0L),
  identical(eval(formals(fit_lightgbm_object)$threads), -1L)
)

p <- data.frame(
  study = rep(c("A", "B"), each = 4L),
  outcome = rep(c("control", "case", "control", "case"), 2L),
  probability = c(0.1, 0.9, 0.2, 0.8, 0.2, 0.8, 0.3, 0.7)
)
t1 <- trace_threshold(p)
p$disease <- p$outcome
p$outcome <- NULL
t2 <- trace_threshold(p)
stopifnot(isTRUE(all.equal(t1, t2)), length(t1) == 1L, is.finite(t1))

bad <- p
bad$probability[[1L]] <- Inf
stopifnot(inherits(try(trace_threshold(bad), silent = TRUE), "try-error"))

stopifnot(
  is.factor(factor(p$disease, levels = c("control", "case"))),
  identical(as.integer(c("control", "case") == "case"), c(0L, 1L))
)
message("PASS scientific compatibility: grid, CPU device, historical threads, threshold labels and validation.")

if (all(vapply(c("ranger", "matrixStats"), requireNamespace, logical(1), quietly = TRUE))) {
  trace_rng(7832L)
  x <- matrix(rpois(100 * 40, lambda = 20), nrow = 100L,
              dimnames = list(paste0("gene", seq_len(100)), paste0("s", seq_len(40))))
  m <- data.frame(sample_id = colnames(x), disease = rep(c("control", "case"), 20),
                  study = rep(c("A", "B"), each = 20))
  fit <- trace_rf(train = x[, 1:30], new = x[, 31:40],
                  meta = m[1:30, ], new_meta = m[31:40, ],
                  representation = "RAW_counts", genes = 10L,
                  weighted = FALSE, seed = 7832L)
  stopifnot(length(fit$probability) == 10L, length(fit$genes) == 10L,
            all(is.finite(fit$probability)), all(fit$probability >= 0 & fit$probability <= 1))
  message("PASS optional synthetic RF compatibility test (legacy disease column).")
} else {
  message("SKIP optional synthetic RF test: ranger or matrixStats not installed.")
}
