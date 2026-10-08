# ======================================================================
# STRICT CCV OF ALL COMPLETE VST + fSVA BENCHMARK CANDIDATES
# ======================================================================

ccv_candidates <- copy(complete_benchmark[Representation_ID == "VST_fSVA"])
ccv_candidates <- ccv_candidates[seq_len(min(CCV_CANDIDATES, .N))]
fwrite(ccv_candidates, file.path(ccv_dir, "TB_strict_CCV_candidates_from_benchmark.csv"))
ccv_rows <- list(); ccv_counter <- 0L
for (candidate_index in seq_len(nrow(ccv_candidates))) {
  candidate <- ccv_candidates[candidate_index]
  for (target_index in seq_along(cohort_names)) {
    target_cohort <- cohort_names[[target_index]]
    source_metadata <- development_metadata[cohort_group != target_cohort]
    target_metadata <- development_metadata[cohort_group == target_cohort]
    source_counts <- development_counts[, source_metadata$sample_id, drop = FALSE]
    target_counts <- development_counts[, target_metadata$sample_id, drop = FALSE]
    result <- tryCatch({
      rep_fit <- fit_apply_vst_fsva(source_counts, target_counts, source_metadata, seed_for(700, candidate_index, target_index))
      fs_fit <- select_top_features(rep_fit$train, rep_fit$new, source_metadata, candidate$FSF_ID[[1L]], TOP_N,
                                    seed_for(710, candidate_index, target_index))
      cba_fit <- apply_cba_training(fs_fit$train, source_metadata, candidate$CBA_ID[[1L]], seed_for(720, candidate_index, target_index))
      source_task <- make_classification_task(cba_fit$matrix, cba_fit$metadata, paste0("CCV_source_", candidate_index, "_", target_index))
      target_task <- make_classification_task(fs_fit$validation, target_metadata, paste0("CCV_target_", candidate_index, "_", target_index))
      learner <- make_tb_learner(candidate$Learner_ID[[1L]], cba_fit$class_weights)
      learner$train(source_task)
      prediction <- learner$predict(target_task)
      metrics <- safe_prediction_metrics(prediction)
      data.table::data.table(
        Candidate_index = candidate_index, Target_cohort = target_cohort,
        Representation = candidate$Representation[[1L]], Representation_ID = candidate$Representation_ID[[1L]],
        FS_method = candidate$FS_method[[1L]], FSF_ID = candidate$FSF_ID[[1L]],
        CBA = candidate$CBA[[1L]], CBA_ID = candidate$CBA_ID[[1L]],
        Learner = candidate$Learner[[1L]], Learner_ID = candidate$Learner_ID[[1L]],
        Source_n = nrow(source_metadata), Target_n = nrow(target_metadata),
        AUROC = metrics$AUROC, AUPRC = metrics$AUPRC, BACC = metrics$BACC, Brier = metrics$Brier,
        Status = "success", Error = NA_character_
      )
    }, error = function(e) data.table::data.table(
      Candidate_index = candidate_index, Target_cohort = target_cohort,
      Representation = candidate$Representation[[1L]], Representation_ID = candidate$Representation_ID[[1L]],
      FS_method = candidate$FS_method[[1L]], FSF_ID = candidate$FSF_ID[[1L]],
      CBA = candidate$CBA[[1L]], CBA_ID = candidate$CBA_ID[[1L]],
      Learner = candidate$Learner[[1L]], Learner_ID = candidate$Learner_ID[[1L]],
      Source_n = nrow(source_metadata), Target_n = nrow(target_metadata),
      AUROC = NA_real_, AUPRC = NA_real_, BACC = NA_real_, Brier = NA_real_,
      Status = "failed", Error = conditionMessage(e)
    ))
    ccv_counter <- ccv_counter + 1L; ccv_rows[[ccv_counter]] <- result
  }
}
ccv_results <- rbindlist(ccv_rows, use.names = TRUE, fill = TRUE)
fwrite(ccv_results, file.path(ccv_dir, "TB_strict_CCV_results.csv"))
ccv_summary <- ccv_results[Status == "success" & is.finite(BACC), .(
  Target_cohorts = .N,
  Mean_CCV_BACC = mean(BACC), SD_CCV_BACC = sd(BACC), Worst_CCV_BACC = min(BACC),
  Mean_CCV_AUROC = mean(AUROC), SD_CCV_AUROC = sd(AUROC), Worst_CCV_AUROC = min(AUROC),
  Mean_CCV_AUPRC = mean(AUPRC), Mean_CCV_Brier = mean(Brier)
), by = .(Candidate_index, Representation, Representation_ID, FS_method, FSF_ID, CBA, CBA_ID, Learner, Learner_ID)]
ccv_summary <- merge(ccv_summary, ccv_candidates[, .(Candidate_index = .I, Mean_BACC, SD_BACC, Worst_BACC, Mean_AUROC, Benchmark_rank)], by = "Candidate_index", all.x = TRUE)
setorder(ccv_summary, -Mean_CCV_BACC, -Worst_CCV_BACC, SD_CCV_BACC, -Mean_BACC, -Mean_CCV_AUROC, SD_BACC)
ccv_summary[, Generalisability_rank := seq_len(.N)]
fwrite(ccv_summary, file.path(ccv_dir, "TB_strict_CCV_summary.csv"))
if (!nrow(ccv_summary)) stop("No strict CCV candidate completed successfully.")
locked_pipeline <- ccv_summary[1L]
fwrite(locked_pipeline, file.path(tuning_dir, "TB_LOCKED_PIPELINE_BEFORE_TUNING.csv"))
print(ccv_summary)

# ======================================================================
# RANDOM HYPERPARAMETER CONFIGURATIONS FOR THE ONE LOCKED PIPELINE
# ======================================================================

sample_tuning_configs <- function(learner_id, n, seed) {
  set.seed(seed)
  if (learner_id == "RBF_SVM") return(data.table::data.table(Config_ID = seq_len(n), cost = exp(runif(n, log(1e-2), log(1e2))), gamma = exp(runif(n, log(1e-5), log(1)))))
  if (learner_id == "Linear_SVM") return(data.table::data.table(Config_ID = seq_len(n), cost = exp(runif(n, log(1e-2), log(1e2)))))
  if (learner_id == "Random_Forest") return(data.table::data.table(Config_ID = seq_len(n), num.trees = sample(200:500, n, replace = TRUE), mtry.ratio = runif(n, 0.10, 0.50), min.node.size = sample(1:8, n, replace = TRUE)))
  stop("No tuning space for learner: ", learner_id)
}

config_to_params <- function(config_row, learner_id) {
  if (learner_id == "RBF_SVM") return(list(cost = config_row$cost[[1L]], gamma = config_row$gamma[[1L]]))
  if (learner_id == "Linear_SVM") return(list(cost = config_row$cost[[1L]]))
  if (learner_id == "Random_Forest") return(list(num.trees = as.integer(config_row$num.trees[[1L]]), mtry.ratio = config_row$mtry.ratio[[1L]], min.node.size = as.integer(config_row$min.node.size[[1L]])))
}

tuning_configs <- sample_tuning_configs(locked_pipeline$Learner_ID[[1L]], TUNING_EVALS, seed_for(800, 1))
fwrite(tuning_configs, file.path(tuning_dir, "TB_tuning_search_space_realised.csv"))

# ======================================================================
# NESTED TUNING — OUTER 5 x INNER 3, OPTIMISE BACC
# preprocessing + FS + CBA are repeated inside every inner training fold.
# ======================================================================

tuning_outer_map <- make_stratified_fold_map(development_metadata, TUNING_OUTER_FOLDS, 1L, seed_for(810, 1))
tuning_archive_rows <- list(); tuning_outer_rows <- list(); archive_counter <- 0L
for (outer_fold in seq_len(TUNING_OUTER_FOLDS)) {
  outer_val_ids <- tuning_outer_map[Fold == outer_fold, sample_id]
  outer_train_ids <- setdiff(development_metadata$sample_id, outer_val_ids)
  outer_train_meta <- development_metadata[match(outer_train_ids, sample_id)]
  outer_val_meta <- development_metadata[match(outer_val_ids, sample_id)]
  inner_map <- make_stratified_fold_map(outer_train_meta, TUNING_INNER_FOLDS, 1L, seed_for(820, outer_fold))
  config_scores <- numeric(nrow(tuning_configs))
  for (config_index in seq_len(nrow(tuning_configs))) {
    inner_bacc <- numeric(TUNING_INNER_FOLDS)
    for (inner_fold in seq_len(TUNING_INNER_FOLDS)) {
      inner_val_ids <- inner_map[Fold == inner_fold, sample_id]
      inner_train_ids <- setdiff(outer_train_meta$sample_id, inner_val_ids)
      inner_train_meta <- outer_train_meta[match(inner_train_ids, sample_id)]
      inner_val_meta <- outer_train_meta[match(inner_val_ids, sample_id)]
      inner_train_counts <- development_counts[, inner_train_meta$sample_id, drop = FALSE]
      inner_val_counts <- development_counts[, inner_val_meta$sample_id, drop = FALSE]
      inner_bacc[[inner_fold]] <- tryCatch({
        rep_fit <- fit_apply_vst_fsva(inner_train_counts, inner_val_counts, inner_train_meta, seed_for(830, outer_fold, config_index, inner_fold))
        fs_fit <- select_top_features(rep_fit$train, rep_fit$new, inner_train_meta, locked_pipeline$FSF_ID[[1L]], TOP_N,
                                      seed_for(840, outer_fold, config_index, inner_fold))
        cba_fit <- apply_cba_training(fs_fit$train, inner_train_meta, locked_pipeline$CBA_ID[[1L]], seed_for(850, outer_fold, config_index, inner_fold))
        tr_task <- make_classification_task(cba_fit$matrix, cba_fit$metadata, paste0("TUNE_INNER_TR_", outer_fold, "_", config_index, "_", inner_fold))
        va_task <- make_classification_task(fs_fit$validation, inner_val_meta, paste0("TUNE_INNER_VA_", outer_fold, "_", config_index, "_", inner_fold))
        learner <- make_tb_learner(locked_pipeline$Learner_ID[[1L]], cba_fit$class_weights, config_to_params(tuning_configs[config_index], locked_pipeline$Learner_ID[[1L]]))
        learner$train(tr_task)
        pred <- learner$predict(va_task)
        as.numeric(pred$score(mlr3::msr("classif.bacc")))
      }, error = function(e) NA_real_)
    }
    config_scores[[config_index]] <- mean(inner_bacc, na.rm = TRUE)
    archive_counter <- archive_counter + 1L
    tuning_archive_rows[[archive_counter]] <- cbind(
      data.table::data.table(Outer_fold = outer_fold, Config_ID = tuning_configs$Config_ID[[config_index]], Mean_inner_BACC = config_scores[[config_index]], SD_inner_BACC = sd(inner_bacc, na.rm = TRUE)),
      tuning_configs[config_index, setdiff(names(tuning_configs), "Config_ID"), with = FALSE]
    )
  }
  best_index <- which.max(config_scores)
  best_config <- tuning_configs[best_index]
  best_inner <- config_scores[[best_index]]
  # outer unbiased evaluation
  outer_train_counts <- development_counts[, outer_train_meta$sample_id, drop = FALSE]
  outer_val_counts <- development_counts[, outer_val_meta$sample_id, drop = FALSE]
  rep_fit <- fit_apply_vst_fsva(outer_train_counts, outer_val_counts, outer_train_meta, seed_for(860, outer_fold))
  fs_fit <- select_top_features(rep_fit$train, rep_fit$new, outer_train_meta, locked_pipeline$FSF_ID[[1L]], TOP_N, seed_for(870, outer_fold))
  cba_fit <- apply_cba_training(fs_fit$train, outer_train_meta, locked_pipeline$CBA_ID[[1L]], seed_for(880, outer_fold))
  tr_task <- make_classification_task(cba_fit$matrix, cba_fit$metadata, paste0("TUNE_OUTER_TR_", outer_fold))
  va_task <- make_classification_task(fs_fit$validation, outer_val_meta, paste0("TUNE_OUTER_VA_", outer_fold))
  learner <- make_tb_learner(locked_pipeline$Learner_ID[[1L]], cba_fit$class_weights, config_to_params(best_config, locked_pipeline$Learner_ID[[1L]]))
  learner$train(tr_task); pred <- learner$predict(va_task); metrics <- safe_prediction_metrics(pred)
  tuning_outer_rows[[outer_fold]] <- cbind(
    data.table::data.table(Outer_fold = outer_fold, Best_Config_ID = best_config$Config_ID[[1L]], Best_inner_BACC = best_inner,
                           Outer_BACC = metrics$BACC, Outer_AUROC = metrics$AUROC, Outer_AUPRC = metrics$AUPRC, Outer_Brier = metrics$Brier,
                           Optimisation_gap = best_inner - metrics$BACC),
    best_config[, setdiff(names(best_config), "Config_ID"), with = FALSE]
  )
}

tuning_archive <- rbindlist(tuning_archive_rows, use.names = TRUE, fill = TRUE)
tuning_outer_results <- rbindlist(tuning_outer_rows, use.names = TRUE, fill = TRUE)
fwrite(tuning_archive, file.path(tuning_dir, "TB_nested_tuning_archive.csv"))
fwrite(tuning_outer_results, file.path(tuning_dir, "TB_nested_tuning_outer_results.csv"))

# Hyperparameter-performance robustness / near-optimality.
tuning_stability <- tuning_archive[, {
  best <- max(Mean_inner_BACC, na.rm = TRUE)
  .(Best_inner_BACC = best, Median_inner_BACC = median(Mean_inner_BACC, na.rm = TRUE),
    Inner_BACC_SD = sd(Mean_inner_BACC, na.rm = TRUE), Inner_BACC_IQR = IQR(Mean_inner_BACC, na.rm = TRUE),
    Inner_BACC_range = max(Mean_inner_BACC, na.rm = TRUE) - min(Mean_inner_BACC, na.rm = TRUE),
    Near_optimal_0.005 = mean(Mean_inner_BACC >= best - 0.005, na.rm = TRUE),
    Near_optimal_0.010 = mean(Mean_inner_BACC >= best - 0.010, na.rm = TRUE),
    Near_optimal_0.020 = mean(Mean_inner_BACC >= best - 0.020, na.rm = TRUE))
}, by = Outer_fold]
fwrite(tuning_stability, file.path(tuning_dir, "TB_hyperparameter_performance_stability.csv"))

# Aggregate the SAME realised configs across outer folds. Lock the robust config by
# mean inner BACC, then near-optimal frequency, then SD.
config_robustness <- tuning_archive[, .(
  Mean_inner_BACC = mean(Mean_inner_BACC, na.rm = TRUE),
  SD_between_outer = sd(Mean_inner_BACC, na.rm = TRUE),
  Worst_outer_inner_BACC = min(Mean_inner_BACC, na.rm = TRUE)
), by = Config_ID]
best_by_outer <- tuning_archive[, .(Best = max(Mean_inner_BACC, na.rm = TRUE)), by = Outer_fold]
tmp <- merge(tuning_archive[, .(Outer_fold, Config_ID, Mean_inner_BACC)], best_by_outer, by = "Outer_fold")
near_freq <- tmp[, .(Near_optimal_frequency_0.01 = mean(Mean_inner_BACC >= Best - 0.01)), by = Config_ID]
config_robustness <- merge(config_robustness, near_freq, by = "Config_ID")
config_robustness <- merge(config_robustness, tuning_configs, by = "Config_ID")
setorder(config_robustness, -Mean_inner_BACC, -Near_optimal_frequency_0.01, -Worst_outer_inner_BACC, SD_between_outer)
config_robustness[, Robust_rank := seq_len(.N)]
robust_config <- config_robustness[1L]
fwrite(config_robustness, file.path(tuning_dir, "TB_hyperparameter_robustness_summary.csv"))
fwrite(robust_config, file.path(tuning_dir, "TB_LOCKED_ROBUST_HYPERPARAMETERS.csv"))
print(tuning_outer_results)
print(config_robustness)

# ======================================================================
# FINAL ICE — ONE LOCKED PIPELINE, ONE LOCKED ROBUST CONFIG, ONE TEST
# ======================================================================

final_rep <- fit_apply_vst_fsva(development_counts, final_ice_counts, development_metadata, seed_for(900, 1))
final_features <- select_top_features(final_rep$train, final_rep$new, development_metadata, locked_pipeline$FSF_ID[[1L]], TOP_N, seed_for(910, 1))
final_cba <- apply_cba_training(final_features$train, development_metadata, locked_pipeline$CBA_ID[[1L]], seed_for(920, 1))
final_train_task <- make_classification_task(final_cba$matrix, final_cba$metadata, "TB_FINAL_TRAIN")
final_test_task <- make_classification_task(final_features$validation, final_ice_metadata, "TB_FINAL_ICE")
final_params <- config_to_params(robust_config, locked_pipeline$Learner_ID[[1L]])
set.seed(seed_for(930, 1))
final_learner <- make_tb_learner(locked_pipeline$Learner_ID[[1L]], final_cba$class_weights, final_params)
final_learner$train(final_train_task)
final_prediction <- final_learner$predict(final_test_task)
final_metrics <- safe_prediction_metrics(final_prediction)
final_truth <- factor(as.character(final_prediction$truth), levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
final_roc <- pROC::roc(final_truth, final_metrics$probability_case, levels = c(NEGATIVE_CLASS, POSITIVE_CLASS), direction = "<", quiet = TRUE)
final_ci <- pROC::ci.auc(final_roc, method = "delong", conf.level = 0.95)

final_result <- data.table::data.table(
  Disease = "TB", Representation = "VST + fSVA", FS_method = locked_pipeline$FS_method[[1L]],
  CBA = locked_pipeline$CBA[[1L]], Model = locked_pipeline$Learner[[1L]],
  Development_n = nrow(development_metadata), Final_ice_n = nrow(final_ice_metadata),
  AUROC = final_metrics$AUROC, AUROC_lower95 = as.numeric(final_ci[[1L]]), AUROC_upper95 = as.numeric(final_ci[[3L]]),
  AUPRC = final_metrics$AUPRC, BACC = final_metrics$BACC, Brier = final_metrics$Brier
)
final_predictions <- data.table::data.table(
  sample_id = final_ice_metadata$sample_id, cohort_group = final_ice_metadata$cohort_group,
  truth = as.character(final_prediction$truth), probability_case = final_metrics$probability_case,
  response = as.character(final_prediction$response)
)
fwrite(final_result, file.path(final_ice_dir, "TB_FINAL_ICE_locked_tuned_pipeline_results.csv"))
fwrite(final_predictions, file.path(final_ice_dir, "TB_FINAL_ICE_locked_tuned_pipeline_predictions.csv"))
fwrite(data.table::data.table(Gene_rank = seq_along(final_features$genes), Gene = final_features$genes), file.path(final_ice_dir, "TB_FINAL_ICE_locked_tuned_pipeline_genes.csv"))

# ======================================================================
# PUBLICATION VALIDATION TRAJECTORY TABLE
# Repeated CV -> strict CCV -> nested CV -> untouched Final Ice
# ======================================================================

locked_benchmark_row <- complete_benchmark[
  Representation_ID == locked_pipeline$Representation_ID[[1L]] &
    FSF_ID == locked_pipeline$FSF_ID[[1L]] &
    CBA_ID == locked_pipeline$CBA_ID[[1L]] &
    Learner_ID == locked_pipeline$Learner_ID[[1L]]
][1L]

if (nrow(locked_benchmark_row) != 1L) {
  stop("Could not recover the locked pipeline from the repeated-CV benchmark summary.")
}

validation_trajectory <- data.table::data.table(
  Disease = "TB",
  Representation = "VST + fSVA",
  Model = locked_pipeline$Learner[[1L]],
  FS_method = locked_pipeline$FS_method[[1L]],
  CBA = locked_pipeline$CBA[[1L]],
  Repeated_CV_AUROC = locked_benchmark_row$Mean_AUROC[[1L]],
  Repeated_CV_SD = locked_benchmark_row$SD_AUROC[[1L]],
  CCV_AUROC = locked_pipeline$Mean_CCV_AUROC[[1L]],
  CCV_SD = locked_pipeline$SD_CCV_AUROC[[1L]],
  Nested_CV_AUROC = mean(tuning_outer_results$Outer_AUROC, na.rm = TRUE),
  Nested_CV_SD = stats::sd(tuning_outer_results$Outer_AUROC, na.rm = TRUE),
  Test_AUROC = final_result$AUROC[[1L]],
  Test_AUROC_lower95 = final_result$AUROC_lower95[[1L]],
  Test_AUROC_upper95 = final_result$AUROC_upper95[[1L]]
)

validation_trajectory[, `:=`(
  `Repeated CV AUROC (SD)` = sprintf(
    "%.3f (%.3f)",
    Repeated_CV_AUROC,
    Repeated_CV_SD
  ),
  `CCV AUROC (SD)` = sprintf(
    "%.3f (%.3f)",
    CCV_AUROC,
    CCV_SD
  ),
  `Nested CV AUROC (SD)` = sprintf(
    "%.3f (%.3f)",
    Nested_CV_AUROC,
    Nested_CV_SD
  ),
  `Test AUROC (95% CI)` = sprintf(
    "%.3f (%.3f--%.3f)",
    Test_AUROC,
    Test_AUROC_lower95,
    Test_AUROC_upper95
  )
)]

fwrite(
  validation_trajectory,
  file.path(final_ice_dir, "TB_validation_trajectory_numeric_and_formatted.csv")
)

fwrite(
  validation_trajectory[, .(
    Disease,
    Model,
    FS_method,
    CBA,
    `Repeated CV AUROC (SD)`,
    `CCV AUROC (SD)`,
    `Nested CV AUROC (SD)`,
    `Test AUROC (95% CI)`
  )],
  file.path(final_ice_dir, "TB_PUBLICATION_validation_trajectory_table.csv")
)

print(
  validation_trajectory[, .(
    Disease,
    Model,
    FS_method,
    CBA,
    `Repeated CV AUROC (SD)`,
    `CCV AUROC (SD)`,
    `Nested CV AUROC (SD)`,
    `Test AUROC (95% CI)`
  )]
)

# Quick frozen-score permutation test.
set.seed(seed_for(940, 1))
null_auc <- numeric(FINAL_PERMUTATIONS)
for (b in seq_len(FINAL_PERMUTATIONS)) {
  perm_truth <- final_predictions$truth
  for (cohort_name in unique(final_predictions$cohort_group)) {
    idx <- which(final_predictions$cohort_group == cohort_name)
    perm_truth[idx] <- sample(perm_truth[idx], length(idx), replace = FALSE)
  }
  null_auc[[b]] <- tryCatch(as.numeric(pROC::auc(pROC::roc(factor(perm_truth, levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)), final_predictions$probability_case,
                                                              levels = c(NEGATIVE_CLASS, POSITIVE_CLASS), direction = "<", quiet = TRUE))), error = function(e) NA_real_)
}
empirical_p <- (1 + sum(null_auc >= final_metrics$AUROC, na.rm = TRUE)) / (FINAL_PERMUTATIONS + 1)
permutation_summary <- data.table::data.table(
  Observed_AUROC = final_metrics$AUROC, Null_mean = mean(null_auc, na.rm = TRUE),
  Null_lower95 = quantile(null_auc, 0.025, na.rm = TRUE), Null_upper95 = quantile(null_auc, 0.975, na.rm = TRUE),
  Permutations = FINAL_PERMUTATIONS, Empirical_p = empirical_p
)
fwrite(permutation_summary, file.path(permutation_dir, "TB_FINAL_ICE_permutation_summary.csv"))
fwrite(data.table::data.table(Null_AUROC = null_auc), file.path(permutation_dir, "TB_FINAL_ICE_permutation_null.csv"))

# ======================================================================
# PUBLICATION MODEL BUNDLE + REPRODUCIBLE PREDICTION HELPER
# ======================================================================

predict_tb_bundle <- function(bundle, new_counts) {
  if (is.null(colnames(new_counts))) stop("new_counts requires sample column names.")
  missing_genes <- setdiff(rownames(bundle$development_counts), rownames(new_counts))
  if (length(missing_genes) > 0L) stop("New counts are missing ", length(missing_genes), " genes required by the bundle.")
  new_counts <- new_counts[rownames(bundle$development_counts), , drop = FALSE]
  rep_fit <- fit_apply_vst_fsva(
    train_counts = bundle$development_counts,
    new_counts = new_counts,
    train_metadata = bundle$development_metadata,
    seed = bundle$preprocessing_seed
  )
  new_matrix <- rep_fit$new[bundle$selected_genes, , drop = FALSE]
  dummy_meta <- data.table::data.table(
    sample_id = colnames(new_matrix),
    outcome = factor(rep(NEGATIVE_CLASS, ncol(new_matrix)), levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
  )
  task <- make_classification_task(new_matrix, dummy_meta, "TB_BUNDLE_PREDICT")
  pred <- bundle$trained_learner$predict(task)
  data.table::data.table(
    sample_id = colnames(new_matrix),
    probability_case = as.numeric(pred$prob[, POSITIVE_CLASS]),
    predicted_class = as.character(pred$response)
  )
}

model_bundle <- list(
  bundle_version = "1.0",
  disease = "TB",
  created = as.character(Sys.time()),
  input_scale = "raw integer gene x sample counts",
  representation = "VST + fSVA",
  fs_method = locked_pipeline$FS_method[[1L]],
  fsf_id = locked_pipeline$FSF_ID[[1L]],
  cba = locked_pipeline$CBA[[1L]],
  cba_id = locked_pipeline$CBA_ID[[1L]],
  model = locked_pipeline$Learner[[1L]],
  learner_id = locked_pipeline$Learner_ID[[1L]],
  hyperparameters = final_params,
  selected_genes = final_features$genes,
  top_n = TOP_N,
  preprocessing_seed = seed_for(900, 1),
  training_seed = seed_for(930, 1),
  development_counts = development_counts,
  development_metadata = development_metadata,
  development_sample_ids = development_metadata$sample_id,
  final_ice_sample_ids = final_ice_metadata$sample_id,
  trained_learner = final_learner,
  class_weights = final_cba$class_weights,
  threshold = 0.5,
  benchmark_selection = locked_pipeline,
  robust_hyperparameter_summary = robust_config,
  nested_tuning_outer_results = tuning_outer_results,
  validation_trajectory = validation_trajectory,
  final_ice_result = final_result,
  predict_function = predict_tb_bundle,
  session_info = utils::capture.output(sessionInfo())
)
saveRDS(model_bundle, file.path(bundle_dir, "TB_FULL_FINAL_TRAINED_MODEL_BUNDLE.rds"), compress = "xz")
writeLines(c(
  "TB FULL-DATA FINAL TRAINED MODEL BUNDLE",
  paste0("Created: ", model_bundle$created),
  paste0("Representation: ", model_bundle$representation),
  paste0("FS method: ", model_bundle$fs_method),
  paste0("CBA: ", model_bundle$cba),
  paste0("Model: ", model_bundle$model),
  paste0("TOP_N: ", model_bundle$top_n),
  paste0("Development n: ", nrow(development_metadata)),
  paste0("Final Ice n: ", nrow(final_ice_metadata)),
  paste0("Final Ice AUROC: ", sprintf("%.3f", final_result$AUROC)),
  "Prediction: bundle <- readRDS(...); predict_tb_bundle(bundle, new_counts)",
  "The bundle stores the frozen development population because fSVA projection is reproduced using development-only information."
), file.path(bundle_dir, "README_TB_MODEL_BUNDLE.txt"))

# ======================================================================
# FIGURES — SAME CORE OUTPUTS AS EQUAL-SIZE + TUNING ROBUSTNESS
# ======================================================================

set1 <- RColorBrewer::brewer.pal(9, "Set1")
representation_order <- c("VST + fSVA")
filter_order <- fsf_labels[fsf_ids]
model_order <- learner_labels[learner_ids]
cba_order <- cba_labels[cba_ids]

# Figure 1: VST + fSVA representation geometry
geometry_dt <- rbindlist(list(
  representation_diagnostics[, .(Representation, Structure = "Cohort", k_star = Cohort_k_star, Mean_silhouette = Cohort_mean_silhouette)],
  representation_diagnostics[, .(Representation, Structure = "Outcome", k_star = Outcome_k_star, Mean_silhouette = Outcome_mean_silhouette)]
))
geometry_dt[, Structure := factor(Structure, levels = c("Cohort", "Outcome"))]
geometry_dt[, point_label := sprintf("%.3f\nk*=%d", Mean_silhouette, as.integer(k_star))]
p_representation_geometry <- ggplot2::ggplot(
  geometry_dt,
  ggplot2::aes(Structure, Mean_silhouette, colour = Structure)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_point(size = 4.5) +
  ggplot2::geom_text(ggplot2::aes(label = point_label), vjust = -0.8, fontface = "bold", show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = stats::setNames(set1[1:2], c("Cohort", "Outcome"))) +
  ggplot2::labs(
    title = "VST + fSVA geometry in the full-data TB development population",
    subtitle = "Outcome separation is contrasted with residual cohort structure",
    x = NULL,
    y = "Mean silhouette width",
    colour = "Structure"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "right")

# Figure 2: Jaccard heatmap
jp <- copy(jaccard_pairwise)
revjp <- jp[, .(Representation, FS_method_1 = FS_method_2, FS_method_2 = FS_method_1, Top_n, Shared_genes, Jaccard_index)]
diagjp <- CJ(Representation = representation_order, FS_method_1 = filter_order, unique = TRUE)[, `:=`(FS_method_2 = FS_method_1, Top_n = TOP_N, Shared_genes = TOP_N, Jaccard_index = 1)]
jp <- unique(rbindlist(list(jp, revjp, diagjp), fill = TRUE), by = c("Representation", "FS_method_1", "FS_method_2"))
jp[, `:=`(FS_method_1 = factor(FS_method_1, levels = filter_order), FS_method_2 = factor(FS_method_2, levels = rev(filter_order)), Representation = factor(Representation, levels = representation_order))]
p_jaccard_heatmap <- ggplot2::ggplot(jp, ggplot2::aes(FS_method_1, FS_method_2, fill = Jaccard_index)) + ggplot2::geom_tile(colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Jaccard_index)), size = 3, fontface = "bold") + ggplot2::facet_wrap(~Representation) +
  ggplot2::scale_fill_gradientn(colours = viridisLite::viridis(256), limits = c(0,1)) + ggplot2::coord_fixed() +
  ggplot2::labs(title = "Feature-selection concordance in full-data TB", x = "FS method", y = "FS method") + ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(face = "bold"))

# Figure 3: repeated-CV benchmark BACC heatmap
bh <- copy(complete_benchmark)
bh[, `:=`(
  FS_method = factor(FS_method, levels = filter_order),
  Learner = factor(Learner, levels = model_order),
  CBA = factor(CBA, levels = cba_order)
)]
p_benchmark_heatmap <- ggplot2::ggplot(
  bh,
  ggplot2::aes(FS_method, CBA, fill = Mean_BACC)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Mean_BACC)), size = 3.0, fontface = "bold") +
  ggplot2::facet_wrap(~ Learner, nrow = 1) +
  ggplot2::scale_fill_gradientn(colours = viridisLite::viridis(256), limits = c(0.5, 1), oob = scales::squish) +
  ggplot2::labs(
    title = "Full-data TB repeated-CV benchmark",
    subtitle = paste0(CV_FOLDS, "-fold CV x ", CV_REPEATS, " repeats; VST + fSVA only"),
    x = "Feature-selection method",
    y = "Class-balancing approach",
    fill = "Mean BACC"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold")
  )

# Figure 4: benchmark mean vs worst BACC robustness
p_pipeline_robustness <- ggplot2::ggplot(complete_benchmark, ggplot2::aes(Mean_BACC, Worst_BACC, colour = CBA, shape = Learner, size = SD_BACC)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") + ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.85) + ggplot2::scale_colour_manual(values = set1[1:3]) + ggplot2::coord_equal(xlim = c(0.4,1), ylim = c(0.4,1)) +
  ggplot2::labs(title = "Benchmark performance and fold robustness", x = "Mean BACC", y = "Worst-fold BACC") + ggplot2::theme_minimal(base_size = 10)

# Figure 5: strict CCV candidate transportability
ccvp <- merge(ccv_results[Status == "success"], ccv_summary[, .(Candidate_index, Generalisability_rank)], by = "Candidate_index")
ccvp[, Pipeline := paste(Learner, FS_method, CBA, sep = " | ")]
p_ccv_transportability <- ggplot2::ggplot(ccvp, ggplot2::aes(BACC, reorder(Pipeline, Generalisability_rank), colour = Target_cohort)) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed") + ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.4)) +
  ggplot2::scale_colour_manual(values = set1[seq_along(cohort_names)]) + ggplot2::labs(title = "Strict cohort transportability of benchmark-selected candidates", x = "Held-out cohort BACC", y = "Pipeline") + ggplot2::theme_minimal(base_size = 10)

# Figure 6: benchmark vs CCV mean BACC
p_cv_vs_ccv <- ggplot2::ggplot(ccv_summary, ggplot2::aes(Mean_BACC, Mean_CCV_BACC, colour = CBA, shape = Learner)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") + ggplot2::geom_point(size = 3.5) + ggplot2::scale_colour_manual(values = set1[1:3]) +
  ggplot2::coord_equal(xlim = c(0.4,1), ylim = c(0.4,1)) + ggplot2::labs(title = "Internal learnability versus strict cohort generalisability", x = "Mean benchmark BACC", y = "Mean strict-CCV BACC") + ggplot2::theme_minimal(base_size = 10)

# Figure 7: MAIN FIGURE — hyperparameter robustness and performance effect
# Panel A shows the distribution of inner-CV BACC for each realised
# configuration across outer folds. Panel B shows how each tuned
# hyperparameter relates to performance. Set1 is used throughout.

config_plot_dt <- merge(
  tuning_archive,
  config_robustness[, .(Config_ID, Robust_rank, Mean_inner_BACC_overall = Mean_inner_BACC)],
  by = "Config_ID",
  all.x = TRUE
)

config_order <- config_robustness[order(Robust_rank), as.character(Config_ID)]
config_plot_dt[, Config_ID_factor := factor(as.character(Config_ID), levels = config_order)]
config_plot_dt[, Outer_fold_factor := factor(Outer_fold)]

outer_fold_colours <- stats::setNames(
  set1[seq_len(TUNING_OUTER_FOLDS)],
  levels(config_plot_dt$Outer_fold_factor)
)

p_hyperparameter_config <- ggplot2::ggplot(
  config_plot_dt,
  ggplot2::aes(Config_ID_factor, Mean_inner_BACC, colour = Outer_fold_factor)
) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = Config_ID_factor),
    outlier.shape = NA,
    colour = "grey35",
    fill = NA,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0), size = 2.4, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = outer_fold_colours, name = "Outer fold") +
  ggplot2::labs(
    title = "A. Robustness across hyperparameter configurations",
    subtitle = "Configurations are ordered by robustness rank",
    x = "Hyperparameter configuration",
    y = "Mean inner-CV BACC"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

parameter_columns <- setdiff(
  names(tuning_configs),
  "Config_ID"
)

parameter_effect_rows <- lapply(
  parameter_columns,
  function(parameter_name) {
    current <- config_plot_dt[, .(
      Outer_fold,
      Outer_fold_factor,
      Config_ID,
      Mean_inner_BACC,
      Parameter_value = as.numeric(get(parameter_name))
    )]
    current[, Parameter := parameter_name]
    current
  }
)

parameter_effect_dt <- rbindlist(parameter_effect_rows, use.names = TRUE, fill = TRUE)
parameter_effect_dt <- parameter_effect_dt[is.finite(Parameter_value) & is.finite(Mean_inner_BACC)]

p_hyperparameter_effect <- ggplot2::ggplot(
  parameter_effect_dt,
  ggplot2::aes(Parameter_value, Mean_inner_BACC, colour = Outer_fold_factor)
) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_point(size = 2.4, alpha = 0.80) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = 1),
    method = "loess",
    se = TRUE,
    colour = "black",
    fill = "grey80",
    linewidth = 0.75,
    inherit.aes = TRUE,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Parameter, scales = "free_x") +
  ggplot2::scale_colour_manual(values = outer_fold_colours, guide = "none") +
  ggplot2::labs(
    title = "B. Hyperparameter values and predictive performance",
    subtitle = "The shape of each response reveals flat versus sensitive optimisation regions",
    x = "Hyperparameter value",
    y = "Mean inner-CV BACC"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold")
  )

p_tuning_robustness <- p_hyperparameter_config / p_hyperparameter_effect +
  patchwork::plot_layout(heights = c(0.9, 1.1), guides = "collect") +
  patchwork::plot_annotation(
    title = "Hyperparameter robustness of the locked TB pipeline",
    subtitle = paste0(
      locked_pipeline$Learner[[1L]], " | ",
      locked_pipeline$FS_method[[1L]], " | ",
      locked_pipeline$CBA[[1L]],
      "; nested CV optimised for BACC"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
  ) &
  ggplot2::theme(legend.position = "right")

# Figure 8: Final Ice ROC + bootstrap band
specificity_grid <- seq(0, 1, length.out = 101L)
ci_se <- pROC::ci.se(final_roc, specificities = specificity_grid, conf.level = 0.95, boot.n = ROC_BOOTSTRAPS, boot.stratified = TRUE)
ci_mat <- as.matrix(ci_se)
roc_ribbon <- data.table::data.table(Specificity = as.numeric(rownames(ci_mat)), Lower = ci_mat[,1], Median = ci_mat[,2], Upper = ci_mat[,3])[, FPR := 1 - Specificity]
roc_curve <- data.table::data.table(FPR = 1 - final_roc$specificities, Sensitivity = final_roc$sensitivities)
p_final_ice_roc <- ggplot2::ggplot() + ggplot2::geom_ribbon(data = roc_ribbon, ggplot2::aes(FPR, ymin = Lower, ymax = Upper), alpha = 0.15) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") + ggplot2::geom_step(data = roc_curve, ggplot2::aes(FPR, Sensitivity), linewidth = 1.2) +
  ggplot2::coord_equal(xlim = c(0,1), ylim = c(0,1)) + ggplot2::labs(title = "Locked tuned TB pipeline on Final Ice", subtitle = paste0("n = ", nrow(final_ice_metadata), "; AUROC = ", sprintf("%.3f", final_result$AUROC), " (95% DeLong CI ", sprintf("%.3f", final_result$AUROC_lower95), "–", sprintf("%.3f", final_result$AUROC_upper95), ")"), x = "False-positive rate", y = "Sensitivity") + ggplot2::theme_minimal(base_size = 11)

# Figure 9: permutation null
p_final_permutation <- ggplot2::ggplot(data.table::data.table(Null_AUROC = null_auc), ggplot2::aes(Null_AUROC)) + ggplot2::geom_histogram(bins = 30, fill = "grey82", colour = "white") +
  ggplot2::geom_vline(xintercept = final_result$AUROC, linewidth = 1.2) + ggplot2::labs(title = "Final-Ice discrimination relative to the permutation null", subtitle = paste0(FINAL_PERMUTATIONS, " cohort-stratified permutations; empirical p = ", format.pval(empirical_p, digits = 3)), x = "Permuted AUROC", y = "Count") + ggplot2::theme_minimal(base_size = 11)

figure_manifest <- list(
  Figure_1_representation_geometry = p_representation_geometry,
  Figure_2_Jaccard_heatmap = p_jaccard_heatmap,
  Figure_3_benchmark_BACC_heatmap = p_benchmark_heatmap,
  Figure_4_pipeline_robustness = p_pipeline_robustness,
  Figure_5_strict_CCV_transportability = p_ccv_transportability,
  Figure_6_benchmark_vs_CCV = p_cv_vs_ccv,
  Figure_7_hyperparameter_robustness = p_tuning_robustness,
  Figure_8_Final_Ice_ROC = p_final_ice_roc,
  Figure_9_Final_Ice_permutation = p_final_permutation
)
for (nm in names(figure_manifest)) {
  if (nm == "Figure_7_hyperparameter_robustness") {
    current_width <- 260
    current_height <- 235
  } else {
    current_width <- 220
    current_height <- 170
  }

  ggplot2::ggsave(
    file.path(figure_dir, paste0(nm, ".png")),
    figure_manifest[[nm]],
    width = current_width,
    height = current_height,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(nm, ".pdf")),
    figure_manifest[[nm]],
    width = current_width,
    height = current_height,
    units = "mm",
    bg = "white"
  )
  saveRDS(figure_manifest[[nm]], file.path(figure_dir, paste0(nm, "_plot.rds")))
}

# ======================================================================
# FINAL AUDITS / README
# ======================================================================

fwrite(benchmark_results[Status != "success"], file.path(audit_dir, "TB_benchmark_failures.csv"))
fwrite(ccv_results[Status != "success"], file.path(audit_dir, "TB_CCV_failures.csv"))
capture.output(sessionInfo(), file = file.path(audit_dir, "TB_sessionInfo.txt"))
writeLines(c(
  "TB FULL-DATA FINAL MODEL EXPERIMENT",
  "",
  paste0("Full aligned n: ", nrow(tb_metadata)),
  paste0("Development n: ", nrow(development_metadata)),
  paste0("Final Ice n: ", nrow(final_ice_metadata)),
  "",
  "Benchmark: 3 representations x 3 FS methods x 3 models x 3 CBAs; 5-fold CV.",
  "CBAs: None, cohort-stratified SMOTE, class weighting.",
  "Benchmark ranking: mean BACC, worst-fold BACC, SD BACC, AUROC, Brier.",
  paste0("Strict CCV candidates: all complete benchmark pipelines (VST + fSVA only; maximum ", CCV_CANDIDATES, ")."),
  paste0("Repeated benchmark CV: ", CV_FOLDS, " folds x ", CV_REPEATS, " repeats."),
  "One pipeline is locked before tuning.",
  paste0("Tuning: ", TUNING_OUTER_FOLDS, " outer folds x ", TUNING_INNER_FOLDS, " inner folds x ", TUNING_EVALS, " shared random configurations; objective = BACC."),
  "Final robust hyperparameters: highest mean inner BACC, near-optimal frequency, worst-case performance, then SD.",
  "Final Ice opened once after pipeline and hyperparameters were locked.",
  paste0("Final model: ", locked_pipeline$Learner[[1L]], " | ", locked_pipeline$FS_method[[1L]], " | ", locked_pipeline$CBA[[1L]]),
  paste0("Final AUROC: ", sprintf("%.3f", final_result$AUROC)),
  paste0("Final BACC: ", sprintf("%.3f", final_result$BACC)),
  paste0("Permutation p: ", format.pval(empirical_p, digits = 3)),
  "",
  "Publication model bundle: 100_MODEL_BUNDLE/TB_FULL_FINAL_TRAINED_MODEL_BUNDLE.rds"
), file.path(output_root, "README_TB_FULL_DATA_FINAL_MODEL.txt"))

cat("\n============================================================\n")
cat("TB FULL-DATA FINAL MODEL COMPLETE\n")
cat("============================================================\n")
print(split_summary)
cat("\nLOCKED PRE-TUNING PIPELINE\n"); print(locked_pipeline)
cat("\nROBUST HYPERPARAMETERS\n"); print(robust_config)
cat("\nVALIDATION TRAJECTORY\n"); print(validation_trajectory[, .(Disease, Model, FS_method, CBA, `Repeated CV AUROC (SD)`, `CCV AUROC (SD)`, `Nested CV AUROC (SD)`, `Test AUROC (95% CI)`)])
cat("\nFINAL ICE RESULT\n"); print(final_result)
cat("\nMODEL BUNDLE\n", file.path(bundle_dir, "TB_FULL_FINAL_TRAINED_MODEL_BUNDLE.rds"), "\n", sep = "")

