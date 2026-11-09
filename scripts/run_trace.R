# Minimal configuration-driven TRACE-Omics single holdout runner.
# This is a supervised binary task demonstration, NOT every dissertation stage.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !identical(args[[1L]], "--config")) {
  stop("Usage: Rscript scripts/run_trace.R --config config/my_task.yml", call. = FALSE)
}
if (!file.exists("R/core.R")) stop("Run this script from the repository root.")
source("R/core.R")
d <- trace_load(args[[2L]])
spec <- d$config$pipeline
if (is.null(spec)) {
  stop("Configuration must contain a pre-specified pipeline (representation, filter, genes, learner).")
}
metadata <- d$metadata
train_i <- which(metadata$partition == "development")
test_i <- which(metadata$partition == "holdout")
# The pipeline must be selected before the untouched holdout is evaluated.
result_dir <- trace_record(d, "holdout_evaluation")
result <- trace_predict_pipeline(
  train = d$counts[, train_i, drop = FALSE],
  new = d$counts[, test_i, drop = FALSE],
  meta = metadata[train_i, , drop = FALSE],
  new_meta = metadata[test_i, , drop = FALSE],
  spec = spec,
  seed = d$config$seed,
  project_root = d$config$project_root,
  config = d$config
)
predictions <- data.frame(sample_id = metadata$sample_id[test_i],
  study = metadata$study[test_i],
  outcome = as.character(metadata$outcome[test_i]),
  probability = result$probability, stringsAsFactors = FALSE)
utils::write.csv(predictions, file.path(result_dir, "holdout_predictions_private.csv"), row.names = FALSE)
utils::write.csv(data.frame(selected_gene = result$genes),
                 file.path(result_dir, "selected_genes.csv"), row.names = FALSE)
metrics <- trace_metrics(predictions$outcome, predictions$probability)
utils::write.csv(data.frame(metric = names(metrics), value = unname(metrics)),
                 file.path(result_dir, "holdout_metrics.csv"), row.names = FALSE)
cat("Completed holdout evaluation; results are PRIVATE at:", result_dir, "\n")
