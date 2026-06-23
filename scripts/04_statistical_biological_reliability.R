# ==============================================================================
# TRACE-Omics
# Script 04: Statistical and biological reliability
#
# Purpose:
#   This script evaluates final selected models under a within-distribution
#   held-out test setting. It trains final champion and comparator models,
#   evaluates discrimination, threshold-dependent performance, calibration,
#   permutation-based feature importance and Reactome pathway enrichment.
#
# Required input structure:
#
#   data/processed/TB/tasks/final_tasks/
#   ├── train_task_TB_VST_SVA_information_gain_SMOTE.rds
#   ├── test_task_TB_VST_SVA_information_gain_SMOTE.rds
#   ├── train_task_TB_VST_SVA_HybridA_SMOTE.rds
#   └── test_task_TB_VST_SVA_HybridA_SMOTE.rds
#
# Optional nested-CV archive structure:
#
#   results/03_benchmarking_ncv_ccv/TB/02_nested_cv/
#   ├── nestedcv_RBF_SVM_task_TB_VST_SVA_information_gain_SMOTE.rds
#   ├── nestedcv_LightGBM_task_TB_VST_SVA_information_gain_SMOTE.rds
#   └── nestedcv_RandomForest_task_TB_VST_SVA_HybridA_SMOTE.rds
#
# Main outputs:
#   - final model production bundles
#   - held-out test prediction table
#   - test metrics with AUROC CIs, Youden thresholds, BACC, Sens and Spec
#   - ROC, PR and calibration plots
#   - DALEX permutation feature-importance tables and plots
#   - consensus biomarker ranking
#   - consensus and model-specific Reactome enrichment outputs
#   - sessionInfo log
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
  library(ggplot2)
  library(pROC)
  library(mlr3verse)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3extralearners)
  library(DALEX)
  library(ingredients)
})

set.seed(7832)

lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

# ==============================================================================
# 1. User settings
# ==============================================================================

disease_id <- "TB"

positive_class <- "case"
negative_class <- "control"

final_task_dir <- file.path(
  "data", "processed", disease_id, "tasks", "final_tasks"
)

nestedcv_dir <- file.path(
  "results", "03_benchmarking_ncv_ccv", disease_id, "02_nested_cv"
)

base_out_dir <- file.path(
  "results", "04_statistical_biological_reliability", disease_id
)

bundle_dir <- file.path(base_out_dir, "production_bundles")
test_dir <- file.path(base_out_dir, "heldout_test")
explain_dir <- file.path(base_out_dir, "feature_importance")
enrich_dir <- file.path(base_out_dir, "reactome_enrichment")
fig_dir <- file.path(base_out_dir, "figures")
summary_dir <- file.path(base_out_dir, "summary")

dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(explain_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Final champion and comparator configurations.
#
# Edit this table only.
# - task_key must match the train/test task filenames after train_task_ and test_task_
# - learner_key must be one of: RBF_SVM, LightGBM, RandomForest
# - nestedcv_rds can be NA if you want to use fallback_params
# ------------------------------------------------------------------------------

jobs <- tribble(
  ~role,        ~label,                                             ~learner_key,   ~task_key,                                  ~nestedcv_rds,
  "Champion",   "RBF SVM | VST + SVA | Info_Gain | SMOTE",           "RBF_SVM",      paste0(disease_id, "_VST_SVA_information_gain_SMOTE"), NA_character_,
  "Comparator", "LightGBM | VST + SVA | Info_Gain | SMOTE",          "LightGBM",     paste0(disease_id, "_VST_SVA_information_gain_SMOTE"), NA_character_,
  "Comparator", "Random Forest | VST + SVA | Info_Gain | SMOTE",     "RandomForest", paste0(disease_id, "_VST_SVA_information_gain_SMOTE"), NA_character_,
  "Comparator", "Random Forest | VST + SVA | Hybrid A | SMOTE",      "RandomForest", paste0(disease_id, "_VST_SVA_HybridA_SMOTE"),          NA_character_
)

# If nestedcv_rds is NA, the script looks automatically in nestedcv_dir using:
# nestedcv_<learner_key>_task_<task_key>.rds

# Fallback fixed parameters are used only when a nested-CV archive cannot be found.
# Replace these with disease-specific values if you are not using nested-CV archives.

fallback_params <- list(
  RBF_SVM = list(
    cost = 32.24861,
    gamma = 0.05145614
  ),
  LightGBM = list(
    learning_rate = 0.1913281,
    num_leaves = 39,
    feature_fraction = 0.8140357
  ),
  RandomForest = list(
    num.trees = 276,
    mtry.ratio = 0.1981858,
    min.node.size = 5
  )
)

# Optional fixed operating values from outer CV.
# These are used for labelling ROC legends and storing the final operating point.
# Set to NA if not available.

outer_operating_tbl <- data.table(
  Configuration = jobs$label,
  J_outer = c(0.969, 0.916, 0.892, 0.889)
)

config_levels <- jobs$label

pal <- c(
  "RBF SVM | VST + SVA | Info_Gain | SMOTE" = "#E41A1C",
  "LightGBM | VST + SVA | Info_Gain | SMOTE" = "#377EB8",
  "Random Forest | VST + SVA | Info_Gain | SMOTE" = "#4DAF4A",
  "Random Forest | VST + SVA | Hybrid A | SMOTE" = "#984EA3"
)

# Ensure every configuration has a colour.
missing_cols <- setdiff(config_levels, names(pal))
if (length(missing_cols) > 0) {
  extra_cols <- RColorBrewer::brewer.pal(max(3, length(missing_cols)), "Set2")[seq_along(missing_cols)]
  pal <- c(pal, setNames(extra_cols, missing_cols))
}

num_permutations_fi <- 50
num_boot <- 500

# ==============================================================================
# 2. Helper functions
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

normalise_params <- function(params, learner_key) {
  if (is.null(params)) {
    return(fallback_params[[learner_key]])
  }

  if (learner_key == "RBF_SVM") {
    list(
      cost = as.numeric(params[["classif.svm.cost"]] %||% params[["cost"]] %||% fallback_params$RBF_SVM$cost),
      gamma = as.numeric(params[["classif.svm.gamma"]] %||% params[["gamma"]] %||% fallback_params$RBF_SVM$gamma)
    )
  } else if (learner_key == "LightGBM") {
    list(
      learning_rate = as.numeric(params[["classif.lightgbm.learning_rate"]] %||% params[["learning_rate"]] %||% fallback_params$LightGBM$learning_rate),
      num_leaves = as.integer(params[["classif.lightgbm.num_leaves"]] %||% params[["num_leaves"]] %||% fallback_params$LightGBM$num_leaves),
      feature_fraction = as.numeric(params[["classif.lightgbm.feature_fraction"]] %||% params[["feature_fraction"]] %||% fallback_params$LightGBM$feature_fraction)
    )
  } else if (learner_key == "RandomForest") {
    list(
      num.trees = as.integer(params[["classif.ranger.num.trees"]] %||% params[["num.trees"]] %||% fallback_params$RandomForest$num.trees),
      mtry.ratio = as.numeric(params[["classif.ranger.mtry.ratio"]] %||% params[["mtry.ratio"]] %||% fallback_params$RandomForest$mtry.ratio),
      min.node.size = as.integer(params[["classif.ranger.min.node.size"]] %||% params[["min.node.size"]] %||% fallback_params$RandomForest$min.node.size)
    )
  } else {
    stop("Unsupported learner_key: ", learner_key, call. = FALSE)
  }
}

find_nestedcv_path <- function(learner_key, task_key, nestedcv_rds) {
  if (!is.na(nestedcv_rds) && file.exists(nestedcv_rds)) {
    return(nestedcv_rds)
  }

  candidate <- file.path(
    nestedcv_dir,
    paste0("nestedcv_", learner_key, "_task_", task_key, ".rds")
  )

  if (file.exists(candidate)) {
    return(candidate)
  }

  NA_character_
}

extract_best_params_from_nestedcv <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }

  rr <- readRDS(path)

  learner_obj <- tryCatch(rr$learners[[1]], error = function(e) NULL)

  if (is.null(learner_obj)) {
    return(NULL)
  }

  archive_data <- tryCatch(
    learner_obj$archive$data,
    error = function(e) NULL
  )

  if (is.null(archive_data)) {
    best_row <- tryCatch(learner_obj$archive$best(), error = function(e) NULL)

    if (!is.null(best_row) && "x_domain" %in% names(best_row)) {
      return(best_row$x_domain[[1]])
    }

    return(NULL)
  }

  score_col <- grep("classif\\.bacc$", names(archive_data), value = TRUE)

  if (length(score_col) == 0) {
    score_col <- grep("bacc", names(archive_data), value = TRUE)
  }

  if (length(score_col) == 0) {
    return(NULL)
  }

  best_idx <- which.max(archive_data[[score_col[1]]])

  if ("x_domain" %in% names(archive_data)) {
    return(archive_data$x_domain[[best_idx]])
  }

  as.list(archive_data[best_idx, ])
}

build_learner <- function(learner_key, params) {
  params <- normalise_params(params, learner_key)

  base_graph <- po("encode", method = "one-hot") %>>%
    po("removeconstants")

  if (learner_key == "RBF_SVM") {
    return(
      as_learner(
        base_graph %>>%
          po("scale") %>>%
          lrn(
            "classif.svm",
            predict_type = "prob",
            kernel = "radial",
            type = "C-classification",
            cost = params$cost,
            gamma = params$gamma
          )
      )
    )
  }

  if (learner_key == "LightGBM") {
    return(
      as_learner(
        base_graph %>>%
          lrn(
            "classif.lightgbm",
            predict_type = "prob",
            learning_rate = params$learning_rate,
            num_leaves = params$num_leaves,
            feature_fraction = params$feature_fraction,
            num_threads = 1,
            verbose = -1
          )
      )
    )
  }

  if (learner_key == "RandomForest") {
    return(
      as_learner(
        base_graph %>>%
          lrn(
            "classif.ranger",
            predict_type = "prob",
            num.trees = params$num.trees,
            mtry.ratio = params$mtry.ratio,
            min.node.size = params$min.node.size,
            num.threads = 1,
            importance = "none"
          )
      )
    )
  }

  stop("Unsupported learner_key: ", learner_key, call. = FALSE)
}

load_clean_task <- function(path) {
  if (!file.exists(path)) {
    stop("Task file not found: ", path, call. = FALSE)
  }

  task <- readRDS(path)$clone(deep = TRUE)

  if ("cohort" %in% task$feature_names) {
    task$select(setdiff(task$feature_names, "cohort"))
  }

  task$col_roles$stratum <- character(0)
  task$col_roles$group <- character(0)

  task
}

predict_case_prob <- function(pred, positive_class = "case") {
  if (!is.null(pred$prob) && positive_class %in% colnames(pred$prob)) {
    return(as.numeric(pred$prob[, positive_class]))
  }

  dt <- as.data.table(pred)

  if ("prob.case" %in% names(dt)) return(as.numeric(dt$prob.case))
  if ("prob_case" %in% names(dt)) return(as.numeric(dt$prob_case))

  stop("Could not extract predicted probability for positive class.", call. = FALSE)
}

calc_metrics_at_youden <- function(df) {
  roc_obj <- pROC::roc(
    response = df$truth,
    predictor = df$prob_case,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- pROC::ci.auc(roc_obj, method = "delong")

  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "specificity", "sensitivity"),
    transpose = FALSE
  )

  threshold <- as.numeric(best$threshold[1])
  sens <- as.numeric(best$sensitivity[1])
  spec <- as.numeric(best$specificity[1])
  bacc <- (sens + spec) / 2
  j <- sens + spec - 1

  data.table(
    Threshold = threshold,
    AUROC = as.numeric(pROC::auc(roc_obj)),
    AUROC_l = as.numeric(auc_ci[1]),
    AUROC_u = as.numeric(auc_ci[3]),
    Sens = sens,
    Spec = spec,
    BACC = bacc,
    J = j
  )
}

bootstrap_metric_ci <- function(df, B = 1000) {
  n <- nrow(df)

  boot_stats <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    d <- df[idx]

    out <- try(calc_metrics_at_youden(d), silent = TRUE)

    if (inherits(out, "try-error")) {
      return(c(Sens = NA, Spec = NA, BACC = NA, J = NA))
    }

    c(
      Sens = out$Sens,
      Spec = out$Spec,
      BACC = out$BACC,
      J = out$J
    )
  })

  data.table(
    Sens_l = quantile(boot_stats["Sens", ], 0.025, na.rm = TRUE),
    Sens_u = quantile(boot_stats["Sens", ], 0.975, na.rm = TRUE),
    Spec_l = quantile(boot_stats["Spec", ], 0.025, na.rm = TRUE),
    Spec_u = quantile(boot_stats["Spec", ], 0.975, na.rm = TRUE),
    BACC_l = quantile(boot_stats["BACC", ], 0.025, na.rm = TRUE),
    BACC_u = quantile(boot_stats["BACC", ], 0.975, na.rm = TRUE),
    J_l = quantile(boot_stats["J", ], 0.025, na.rm = TRUE),
    J_u = quantile(boot_stats["J", ], 0.975, na.rm = TRUE)
  )
}

# ==============================================================================
# 3. Train final learners and save production bundles
# ==============================================================================

bundle_log <- list()
prediction_rows <- list()

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]

  message("\n====================================================")
  message("Training final model: ", job$role, " | ", job$label)

  train_task_path <- file.path(
    final_task_dir,
    paste0("train_task_", job$task_key, ".rds")
  )

  test_task_path <- file.path(
    final_task_dir,
    paste0("test_task_", job$task_key, ".rds")
  )

  nested_path <- find_nestedcv_path(
    learner_key = job$learner_key,
    task_key = job$task_key,
    nestedcv_rds = job$nestedcv_rds
  )

  raw_params <- extract_best_params_from_nestedcv(nested_path)
  final_params <- normalise_params(raw_params, job$learner_key)

  task_train <- load_clean_task(train_task_path)
  task_test <- load_clean_task(test_task_path)

  learner <- build_learner(job$learner_key, final_params)

  set.seed(7832)
  learner$train(task_train)

  common_features <- intersect(task_train$feature_names, task_test$feature_names)

  if (length(common_features) == 0) {
    stop("No common features between train and test task for: ", job$label, call. = FALSE)
  }

  task_test$select(common_features)

  pred <- learner$predict(task_test)

  pred_dt <- data.table(
    truth = factor(as.character(pred$truth), levels = c(negative_class, positive_class)),
    response = factor(as.character(pred$response), levels = c(negative_class, positive_class)),
    prob_case = predict_case_prob(pred, positive_class),
    Configuration = job$label,
    Role = job$role,
    Learner_key = job$learner_key,
    task_key = job$task_key
  )

  prediction_rows[[job$label]] <- pred_dt

  bundle <- list(
    learner = learner,
    feature_names = task_train$feature_names,
    target_name = task_train$target_names,
    class_levels = task_train$class_names,
    params = final_params,
    meta = list(
      disease = disease_id,
      role = job$role,
      label = job$label,
      learner_key = job$learner_key,
      task_key = job$task_key,
      train_task_path = train_task_path,
      test_task_path = test_task_path,
      nestedcv_path = nested_path,
      saved_at = Sys.time()
    )
  )

  bundle_path <- file.path(
    bundle_dir,
    paste0("bundle_", safe_name(job$label), ".rds")
  )

  saveRDS(bundle, bundle_path)

  bundle_log[[job$label]] <- data.table(
    Role = job$role,
    Configuration = job$label,
    Learner_key = job$learner_key,
    task_key = job$task_key,
    Bundle = bundle_path,
    NestedCV = nested_path
  )
}

bundle_log_tbl <- rbindlist(bundle_log, fill = TRUE)
flat_dt_test <- rbindlist(prediction_rows, fill = TRUE)

flat_dt_test[, Configuration := factor(Configuration, levels = config_levels)]

fwrite(
  bundle_log_tbl,
  file.path(bundle_dir, paste0(disease_id, "_bundle_log.csv"))
)

fwrite(
  flat_dt_test,
  file.path(test_dir, paste0(disease_id, "_final_test_prediction_level.csv"))
)

saveRDS(
  flat_dt_test,
  file.path(test_dir, paste0(disease_id, "_final_test_prediction_level.rds"))
)

# ==============================================================================
# 4. Held-out test metrics with Youden threshold optimisation
# ==============================================================================

set.seed(7832)

test_metrics_raw <- rbindlist(
  lapply(
    split(flat_dt_test, flat_dt_test$Configuration),
    function(df) {
      config <- unique(as.character(df$Configuration))
      point <- calc_metrics_at_youden(as.data.table(df))
      ci <- bootstrap_metric_ci(as.data.table(df), B = 1000)

      cbind(
        data.table(Configuration = config),
        point,
        ci
      )
    }
  ),
  fill = TRUE
)

test_metrics_raw <- merge(
  test_metrics_raw,
  outer_operating_tbl,
  by = "Configuration",
  all.x = TRUE
)

test_metrics_publish <- test_metrics_raw[, .(
  Configuration,
  Threshold = sprintf("%.3f", Threshold),
  AUROC = sprintf("%.3f [%.3f--%.3f]", AUROC, AUROC_l, AUROC_u),
  BACC = sprintf("%.3f [%.3f--%.3f]", BACC, BACC_l, BACC_u),
  Sensitivity = sprintf("%.3f [%.3f--%.3f]", Sens, Sens_l, Sens_u),
  Specificity = sprintf("%.3f [%.3f--%.3f]", Spec, Spec_l, Spec_u),
  `J-index` = sprintf("%.3f [%.3f--%.3f]", J, J_l, J_u),
  `Outer-CV J` = ifelse(is.na(J_outer), NA_character_, sprintf("%.3f", J_outer))
)]

fwrite(
  test_metrics_raw,
  file.path(test_dir, paste0(disease_id, "_final_test_metrics_raw.csv"))
)

fwrite(
  test_metrics_publish,
  file.path(test_dir, paste0(disease_id, "_final_test_metrics_publish.csv"))
)

print(test_metrics_publish)

# ==============================================================================
# 5. Final test ROC plot
# ==============================================================================

fpr_grid <- seq(0, 1, by = 0.01)

roc_results <- lapply(
  split(flat_dt_test, flat_dt_test$Configuration),
  function(df) {
    config_name <- as.character(unique(df$Configuration))

    roc_obj <- pROC::roc(
      response = df$truth,
      predictor = df$prob_case,
      levels = c(negative_class, positive_class),
      direction = "<",
      quiet = TRUE
    )

    raw_roc <- as.data.table(
      pROC::coords(
        roc_obj,
        x = "all",
        ret = c("specificity", "sensitivity"),
        transpose = FALSE
      )
    )

    raw_roc[, fpr := 1 - as.numeric(specificity)]
    raw_roc[, sensitivity := as.numeric(sensitivity)]

    raw_roc <- rbindlist(list(
      data.table(fpr = 0, sensitivity = 0),
      raw_roc[, .(fpr, sensitivity)],
      data.table(fpr = 1, sensitivity = 1)
    ), fill = TRUE)

    setorder(raw_roc, fpr, sensitivity)
    raw_roc <- unique(raw_roc, by = "fpr")
    raw_roc[, sensitivity := cummax(sensitivity)]

    interp_sens <- approx(
      x = raw_roc$fpr,
      y = raw_roc$sensitivity,
      xout = fpr_grid,
      method = "linear",
      ties = "ordered",
      rule = 2
    )$y

    roc_df <- data.table(
      fpr = fpr_grid,
      sensitivity = pmin(pmax(interp_sens, 0), 1),
      Configuration = config_name
    )

    ci_obj <- try(
      pROC::ci.se(
        roc_obj,
        specificities = 1 - fpr_grid,
        boot.n = num_boot
      ),
      silent = TRUE
    )

    if (!inherits(ci_obj, "try-error")) {
      ci_tmp <- as.data.table(as.data.frame(ci_obj))
      setnames(ci_tmp, c("sens_lower", "sensitivity_mid", "sens_upper"))

      ci_df <- data.table(
        fpr = fpr_grid,
        sens_lower = pmax(as.numeric(ci_tmp$sens_lower), 0),
        sens_upper = pmin(as.numeric(ci_tmp$sens_upper), 1),
        Configuration = config_name
      )
    } else {
      ci_df <- copy(roc_df)
      ci_df[, `:=`(
        sens_lower = sensitivity,
        sens_upper = sensitivity
      )]
    }

    auc_ci <- pROC::ci.auc(roc_obj, method = "delong")

    auc_df <- data.table(
      Configuration = config_name,
      auc = as.numeric(pROC::auc(roc_obj)),
      low = as.numeric(auc_ci[1]),
      high = as.numeric(auc_ci[3])
    )

    list(
      roc_df = roc_df,
      ci_df = ci_df,
      auc_df = auc_df
    )
  }
)

roc_df_all <- rbindlist(lapply(roc_results, `[[`, "roc_df"), fill = TRUE)
ci_df_all <- rbindlist(lapply(roc_results, `[[`, "ci_df"), fill = TRUE)
auc_df_all <- rbindlist(lapply(roc_results, `[[`, "auc_df"), fill = TRUE)

config_levels_roc <- auc_df_all[order(-auc), as.character(Configuration)]

roc_df_all[, Configuration := factor(Configuration, levels = config_levels_roc)]
ci_df_all[, Configuration := factor(Configuration, levels = config_levels_roc)]
auc_df_all[, Configuration := factor(Configuration, levels = config_levels_roc)]

auc_df_all <- merge(
  auc_df_all,
  outer_operating_tbl,
  by = "Configuration",
  all.x = TRUE
)

auc_df_all[, legend_label := ifelse(
  is.na(J_outer),
  sprintf("%s (AUROC %.3f [%.3f--%.3f])", Configuration, auc, low, high),
  sprintf("%s (AUROC %.3f [%.3f--%.3f]; J = %.3f)", Configuration, auc, low, high, J_outer)
)]

legend_labels_roc <- setNames(
  auc_df_all$legend_label,
  as.character(auc_df_all$Configuration)
)

p_roc_test <- ggplot() +
  geom_ribbon(
    data = ci_df_all,
    aes(x = fpr, ymin = sens_lower, ymax = sens_upper, fill = Configuration),
    alpha = 0.14,
    colour = NA
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.7
  ) +
  geom_line(
    data = roc_df_all,
    aes(x = fpr, y = sensitivity, colour = Configuration),
    linewidth = 1.15,
    lineend = "round"
  ) +
  scale_colour_manual(
    values = pal,
    breaks = config_levels_roc,
    labels = legend_labels_roc[config_levels_roc],
    name = "Configuration"
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.035, 0.045))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    title = paste0(disease_id, " final test ROC comparison"),
    subtitle = "Champion versus scientific comparators",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 8),
    axis.title = element_text(colour = "black", face = "bold", size = 10.5),
    legend.title = element_text(face = "bold", colour = "black", size = 9),
    legend.text = element_text(colour = "black", size = 7.2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.margin = margin(8, 22, 8, 22)
  ) +
  guides(colour = guide_legend(ncol = 1))

print(p_roc_test)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_ROC_champion_vs_comparators.png")),
  plot = p_roc_test,
  width = 240,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_ROC_champion_vs_comparators.pdf")),
  plot = p_roc_test,
  width = 240,
  height = 140,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 6. Final test PR plot
# ==============================================================================

fig_test_dt <- as_tibble(flat_dt_test) %>%
  mutate(
    y = ifelse(as.character(truth) == positive_class, 1, 0),
    p = prob_case,
    Configuration = factor(as.character(Configuration), levels = config_levels)
  ) %>%
  filter(!is.na(p), !is.na(Configuration)) %>%
  group_by(Configuration) %>%
  filter(n_distinct(y) == 2) %>%
  ungroup()

auprc_once <- function(y, p) {
  ord <- order(-p)
  y <- y[ord]

  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)

  recall <- tp / sum(y == 1)
  precision <- tp / pmax(tp + fp, 1)

  recall <- c(0, recall)
  precision <- c(1, precision)

  sum(diff(recall) * precision[-length(precision)], na.rm = TRUE)
}

pr_points_one <- function(df) {
  ord <- order(-df$p)
  y <- df$y[ord]

  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)

  recall <- tp / sum(y == 1)
  precision <- tp / pmax(tp + fp, 1)

  tibble(
    recall = c(0, recall),
    precision = c(precision[1], precision)
  ) %>%
    arrange(recall, desc(precision)) %>%
    distinct(recall, .keep_all = TRUE)
}

boot_pr_one <- function(df, grid = seq(0, 1, by = 0.01), B = num_boot) {
  n <- nrow(df)

  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    d_b <- df[idx, , drop = FALSE]
    pr_b <- pr_points_one(d_b)

    approx(
      x = pr_b$recall,
      y = pr_b$precision,
      xout = grid,
      rule = 2,
      ties = "ordered"
    )$y
  })

  boot_auprc <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    auprc_once(df$y[idx], df$p[idx])
  })

  list(
    ribbon = tibble(
      recall = grid,
      lo = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
      hi = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    ),
    auprc = auprc_once(df$y, df$p),
    low = quantile(boot_auprc, 0.025, na.rm = TRUE),
    high = quantile(boot_auprc, 0.975, na.rm = TRUE),
    prevalence = mean(df$y)
  )
}

split_dt <- split(fig_test_dt, fig_test_dt$Configuration)

pr_points <- list()
pr_ribbon <- list()
pr_table <- list()

for (nm in names(split_dt)) {
  df <- split_dt[[nm]]
  boot <- boot_pr_one(df, B = num_boot)

  pr_points[[nm]] <- pr_points_one(df) %>%
    mutate(Configuration = nm)

  pr_ribbon[[nm]] <- boot$ribbon %>%
    mutate(Configuration = nm)

  pr_table[[nm]] <- tibble(
    Configuration = nm,
    AUPRC = as.numeric(boot$auprc),
    low = as.numeric(boot$low),
    high = as.numeric(boot$high),
    prevalence = as.numeric(boot$prevalence),
    legend_label = sprintf(
      "%s (AUPRC %.3f [%.3f--%.3f], Prev: %.0f%%)",
      nm,
      boot$auprc,
      boot$low,
      boot$high,
      boot$prevalence * 100
    )
  )
}

pr_points <- bind_rows(pr_points)
pr_ribbon <- bind_rows(pr_ribbon)
pr_table <- bind_rows(pr_table)

config_levels_pr <- pr_table %>%
  arrange(desc(AUPRC)) %>%
  pull(Configuration) %>%
  as.character()

pr_points$Configuration <- factor(pr_points$Configuration, levels = config_levels_pr)
pr_ribbon$Configuration <- factor(pr_ribbon$Configuration, levels = config_levels_pr)
pr_table$Configuration <- factor(pr_table$Configuration, levels = config_levels_pr)

fwrite(
  as.data.table(pr_table),
  file.path(test_dir, paste0(disease_id, "_final_test_PR_AUPRC_bootstrap95CI_table.csv"))
)

legend_labels_pr <- setNames(
  pr_table$legend_label,
  as.character(pr_table$Configuration)
)

p_pr_test <- ggplot() +
  geom_hline(
    data = pr_table,
    aes(yintercept = prevalence, colour = Configuration),
    linetype = "dashed",
    linewidth = 0.65,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = pr_ribbon,
    aes(x = recall, ymin = lo, ymax = hi, fill = Configuration),
    alpha = 0.14,
    colour = NA
  ) +
  geom_line(
    data = pr_points,
    aes(x = recall, y = precision, colour = Configuration),
    linewidth = 1.15,
    lineend = "round",
    linejoin = "round"
  ) +
  scale_colour_manual(
    values = pal,
    breaks = config_levels_pr,
    labels = legend_labels_pr[config_levels_pr],
    name = "Configuration (AUPRC [95% CI], prevalence)"
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.035, 0.045))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(
    title = paste0(disease_id, " final test precision-recall comparison"),
    subtitle = "Champion versus scientific comparators; dashed lines indicate case prevalence",
    x = "Recall (Sensitivity)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 8),
    axis.title = element_text(colour = "black", face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5),
    plot.subtitle = element_text(colour = "black", size = 10.5, hjust = 0.5),
    legend.title = element_text(face = "bold", colour = "black", size = 9),
    legend.text = element_text(colour = "black", size = 7.2),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 22, 8, 22)
  ) +
  guides(colour = guide_legend(ncol = 1))

print(p_pr_test)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_PR_champion_vs_comparators.png")),
  plot = p_pr_test,
  width = 240,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_PR_champion_vs_comparators.pdf")),
  plot = p_pr_test,
  width = 240,
  height = 140,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 7. Final test calibration plot
# ==============================================================================

loess_span <- 0.75
belt_boot <- 500
belt_grid <- seq(0, 1, by = 0.01)

cal_test_dt <- as_tibble(flat_dt_test) %>%
  mutate(
    y = ifelse(as.character(truth) == positive_class, 1, 0),
    p = prob_case,
    model = factor(as.character(Configuration), levels = config_levels)
  ) %>%
  filter(!is.na(p), !is.na(model)) %>%
  group_by(model) %>%
  filter(n_distinct(y) == 2) %>%
  ungroup()

bin_summaries_all <- cal_test_dt %>%
  group_by(model) %>%
  mutate(bin = ntile(p, 10)) %>%
  group_by(model, bin) %>%
  summarise(
    n = n(),
    conf = mean(p),
    acc = mean(y),
    acc_lo = binom.test(sum(y), n)$conf.int[1],
    acc_hi = binom.test(sum(y), n)$conf.int[2],
    .groups = "drop"
  )

loess_lines_all <- cal_test_dt %>%
  group_by(model) %>%
  group_modify(~ {
    fit <- loess(y ~ p, span = loess_span, data = .x)

    tibble(
      p = belt_grid,
      m = pmin(
        pmax(
          as.numeric(predict(fit, newdata = data.frame(p = belt_grid))),
          0
        ),
        1
      )
    )
  }) %>%
  ungroup()

belt_df_all <- cal_test_dt %>%
  group_by(model) %>%
  group_modify(~ {
    n <- nrow(.x)

    mats <- matrix(
      NA_real_,
      nrow = length(belt_grid),
      ncol = belt_boot
    )

    for (b in seq_len(belt_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      d_b <- .x[idx, ]

      fit <- try(
        loess(y ~ p, span = loess_span, data = d_b),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        mats[, b] <- as.numeric(
          predict(
            fit,
            newdata = data.frame(p = belt_grid)
          )
        )
      }
    }

    qs <- t(
      apply(
        mats,
        1,
        function(v) {
          quantile(
            v[is.finite(v)],
            probs = c(0.025, 0.975),
            na.rm = TRUE
          )
        }
      )
    )

    tibble(
      p = belt_grid,
      lo = pmax(pmin(qs[, 1], 1), 0),
      hi = pmax(pmin(qs[, 2], 1), 0)
    )
  }) %>%
  ungroup()

ece_metrics <- bin_summaries_all %>%
  group_by(model) %>%
  summarise(
    ECE = sum((n / sum(n)) * abs(acc - conf)),
    .groups = "drop"
  )

brier_metrics <- cal_test_dt %>%
  group_by(model) %>%
  summarise(
    Brier = mean((p - y)^2),
    .groups = "drop"
  )

metric_leg_all <- brier_metrics %>%
  left_join(ece_metrics, by = "model") %>%
  mutate(
    leg_label = sprintf(
      "%s (Brier %.3f, ECE %.3f)",
      model,
      Brier,
      ECE
    )
  )

plot_df_all <- bin_summaries_all %>%
  left_join(metric_leg_all, by = "model")

loess_lines_all <- loess_lines_all %>%
  left_join(metric_leg_all, by = "model")

belt_df_all <- belt_df_all %>%
  left_join(metric_leg_all, by = "model")

legend_labels_cal <- setNames(
  metric_leg_all$leg_label,
  metric_leg_all$model
)

fwrite(
  as.data.table(metric_leg_all),
  file.path(test_dir, paste0(disease_id, "_final_test_calibration_metrics.csv"))
)

p_cal_test <- ggplot() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.8
  ) +
  geom_ribbon(
    data = belt_df_all,
    aes(x = p, ymin = lo, ymax = hi, fill = model),
    alpha = 0.10
  ) +
  geom_line(
    data = loess_lines_all,
    aes(x = p, y = m, colour = model),
    linewidth = 1.3
  ) +
  geom_point(
    data = plot_df_all,
    aes(x = conf, y = acc, colour = model),
    size = 2.3,
    alpha = 0.9
  ) +
  geom_errorbar(
    data = plot_df_all,
    aes(x = conf, ymin = acc_lo, ymax = acc_hi, colour = model),
    width = 0.015,
    linewidth = 0.5,
    alpha = 0.8
  ) +
  scale_x_continuous(
    "Predicted probability of case",
    labels = label_percent(),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    "Observed fraction of cases",
    labels = label_percent(),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_colour_manual(
    values = pal,
    labels = legend_labels_cal,
    name = "Configuration (Brier, ECE)"
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  coord_fixed() +
  labs(
    title = paste0(disease_id, " final test calibration comparison"),
    subtitle = "Champion versus scientific comparators"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    panel.grid.major = element_line(color = "grey92"),
    panel.border = element_rect(linewidth = 1.0, fill = NA),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_cal_test)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_calibration_champion_vs_comparators.png")),
  plot = p_cal_test,
  width = 240,
  height = 140,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_final_test_calibration_champion_vs_comparators.pdf")),
  plot = p_cal_test,
  width = 240,
  height = 140,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 8. DALEX permutation feature importance
# ==============================================================================

pred_fun <- function(model, newdata) {
  nd <- as.data.frame(newdata)
  nd$.target <- factor(negative_class, levels = c(negative_class, positive_class))

  task_new <- TaskClassif$new(
    id = "explain_task",
    backend = nd,
    target = ".target",
    positive = positive_class
  )

  pred <- model$predict(task_new)
  as.numeric(pred$prob[, positive_class])
}

loss_one_minus_auc <- function(observed, predicted, ...) {
  observed <- factor(
    ifelse(observed == 1, positive_class, negative_class),
    levels = c(negative_class, positive_class)
  )

  1 - as.numeric(
    pROC::auc(
      observed,
      predicted,
      levels = c(negative_class, positive_class),
      direction = "<",
      quiet = TRUE
    )
  )
}

all_fi <- list()

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]

  message("\n====================================================")
  message("Permutation importance: ", job$role, " | ", job$label)

  bundle_path <- file.path(
    bundle_dir,
    paste0("bundle_", safe_name(job$label), ".rds")
  )

  if (!file.exists(bundle_path)) {
    warning("Missing bundle: ", bundle_path)
    next
  }

  bundle <- readRDS(bundle_path)
  learner <- bundle$learner

  test_task_path <- file.path(
    final_task_dir,
    paste0("test_task_", job$task_key, ".rds")
  )

  task_test <- load_clean_task(test_task_path)

  common_features <- intersect(bundle$feature_names, task_test$feature_names)

  if (length(common_features) == 0) {
    stop("No common features found for: ", job$label, call. = FALSE)
  }

  task_test$select(common_features)

  x_test <- as.data.frame(task_test$data(cols = common_features))
  y_test <- factor(as.character(task_test$truth()), levels = c(negative_class, positive_class))
  y_num <- as.numeric(y_test == positive_class)

  expl <- DALEX::explain(
    model = learner,
    data = x_test,
    y = y_num,
    predict_function = pred_fun,
    label = job$label,
    verbose = FALSE
  )

  set.seed(7832)

  fi_auc <- ingredients::feature_importance(
    expl,
    loss_function = loss_one_minus_auc,
    type = "difference",
    B = num_permutations_fi
  )

  fi_tbl <- as.data.table(fi_auc)

  fi_tbl[, `:=`(
    Role = job$role,
    Configuration = job$label
  )]

  all_fi[[job$label]] <- fi_tbl

  fwrite(
    fi_tbl,
    file.path(explain_dir, paste0("FI_raw_", safe_name(job$label), ".csv"))
  )
}

fi_all <- rbindlist(all_fi, fill = TRUE)

fwrite(
  fi_all,
  file.path(explain_dir, paste0(disease_id, "_FI_raw_champion_vs_comparators.csv"))
)

fi_clean <- fi_all[
  !variable %in% c("_baseline_", "_full_model_", "_full_model_ ")
]

fi_summary <- fi_clean[
  ,
  .(
    mean_loss = mean(dropout_loss, na.rm = TRUE),
    sd_loss = sd(dropout_loss, na.rm = TRUE)
  ),
  by = .(Role, Configuration, variable)
][order(Configuration, -mean_loss)]

fwrite(
  fi_summary,
  file.path(explain_dir, paste0(disease_id, "_FI_summary_champion_vs_comparators.csv"))
)

gene_rank <- fi_summary[
  ,
  .(
    mean_importance = mean(mean_loss, na.rm = TRUE),
    sd_importance = sd(mean_loss, na.rm = TRUE),
    max_importance = max(mean_loss, na.rm = TRUE),
    n_models = uniqueN(Configuration),
    models_present = paste(unique(Configuration), collapse = "; ")
  ),
  by = variable
][order(-mean_importance)]

gene_rank[, rank := seq_len(.N)]

fwrite(
  gene_rank,
  file.path(explain_dir, paste0(disease_id, "_FI_global_gene_ranking_all_models.csv"))
)

core_biomarkers <- gene_rank[
  n_models >= 2
][order(-mean_importance)]

fwrite(
  core_biomarkers,
  file.path(explain_dir, paste0(disease_id, "_FI_core_biomarker_candidates_consensus.csv"))
)

top_n <- 30

fi_top <- fi_summary[
  ,
  head(.SD[order(-mean_loss)], top_n),
  by = .(Role, Configuration)
]

fi_top[, variable_ordered := forcats::fct_reorder(variable, mean_loss)]
fi_top[, Short_Config := as.character(Configuration)]

p_fi_all <- ggplot(
  fi_top,
  aes(x = variable_ordered, y = mean_loss)
) +
  geom_col(fill = "#377EB8", width = 0.7) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_loss - sd_loss, 0),
      ymax = mean_loss + sd_loss
    ),
    width = 0.2,
    linewidth = 0.3,
    colour = "black",
    alpha = 0.7
  ) +
  coord_flip() +
  facet_wrap(~ Short_Config, scales = "free_y", ncol = 2) +
  scale_y_continuous(
    limits = c(0, max(fi_top$mean_loss, na.rm = TRUE) * 1.08),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Permutation-based feature importance: champion versus comparators",
    subtitle = paste0("Decrease in AUROC following feature permutation on the held-out ", disease_id, " test set"),
    x = "Gene",
    y = expression("Permutation importance (" * Delta * "AUROC after feature permutation)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 8, colour = "black"),
    strip.text = element_text(face = "bold", size = 8.5, colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.7),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_fi_all)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_FI_champion_vs_comparators_top30.png")),
  plot = p_fi_all,
  width = 260,
  height = 210,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, paste0(disease_id, "_FI_champion_vs_comparators_top30.pdf")),
  plot = p_fi_all,
  width = 260,
  height = 210,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

top_consensus_n <- 15

gene_rank_top <- gene_rank[
  n_models >= 2
][order(-mean_importance)][1:min(.N, top_consensus_n)]

if (nrow(gene_rank_top) > 0) {
  gene_rank_top[, variable_ordered := forcats::fct_reorder(variable, mean_importance)]

  p_gene_rank <- ggplot(
    gene_rank_top,
    aes(x = variable_ordered, y = mean_importance)
  ) +
    geom_col(fill = "#377EB8", width = 0.65) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_importance - sd_importance, 0),
        ymax = mean_importance + sd_importance
      ),
      width = 0.2,
      linewidth = 0.25,
      colour = "black",
      alpha = 0.45
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = paste0("Consensus ", disease_id, " biomarker candidates ranked by permutation importance"),
      subtitle = "Genes appearing among important features in at least two final model configurations",
      x = "Gene",
      y = expression("Mean permutation importance (" * Delta * "AUROC)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text.y = element_text(face = "italic", size = 9, colour = "black"),
      axis.text.x = element_text(size = 8, colour = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.7),
      plot.margin = margin(10, 15, 10, 20)
    )

  print(p_gene_rank)

  ggsave(
    filename = file.path(fig_dir, paste0(disease_id, "_FI_consensus_biomarker_candidates_top15.png")),
    plot = p_gene_rank,
    width = 180,
    height = 150,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(fig_dir, paste0(disease_id, "_FI_consensus_biomarker_candidates_top15.pdf")),
    plot = p_gene_rank,
    width = 180,
    height = 150,
    units = "mm",
    device = cairo_pdf,
    bg = "white"
  )
}

# ==============================================================================
# 9. ReactomePA pathway enrichment
# ==============================================================================

required_bioc <- c("org.Hs.eg.db", "AnnotationDbi", "ReactomePA", "clusterProfiler")

missing_bioc <- required_bioc[
  !vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  warning(
    "Skipping Reactome enrichment because these packages are not installed: ",
    paste(missing_bioc, collapse = ", ")
  )
} else {

  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(clusterProfiler)
    library(ReactomePA)
  })

  # ---------------------------------------------------------------------------
  # 9.1 Consensus Reactome enrichment
  # ---------------------------------------------------------------------------

  consensus_genes <- gene_rank[
    n_models >= 2
  ][
    order(-n_models, -mean_importance)
  ][
    1:min(50, .N),
    variable
  ]

  if (length(consensus_genes) < 5) {
    warning("Fewer than 5 consensus genes found. Using top 100 genes by mean importance.")

    consensus_genes <- gene_rank[
      order(-mean_importance)
    ][
      1:min(100, .N),
      variable
    ]
  }

  consensus_genes <- unique(consensus_genes)

  fwrite(
    data.table(SYMBOL = consensus_genes),
    file.path(enrich_dir, paste0(disease_id, "_consensus_genes_used_for_reactome.csv"))
  )

  map_consensus <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = consensus_genes,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    as.data.table() %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, ENTREZID)

  fwrite(
    map_consensus,
    file.path(enrich_dir, paste0(disease_id, "_consensus_gene_symbol_to_entrez.csv"))
  )

  entrez_consensus <- unique(map_consensus$ENTREZID)

  if (length(entrez_consensus) >= 5) {
    reactome_consensus <- ReactomePA::enrichPathway(
      gene = entrez_consensus,
      organism = "human",
      pvalueCutoff = 0.10,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.20,
      readable = TRUE
    )

    if (!is.null(reactome_consensus) && nrow(as.data.frame(reactome_consensus)) > 0) {
      reactome_dt <- as.data.table(as.data.frame(reactome_consensus))

      fwrite(
        reactome_dt,
        file.path(enrich_dir, paste0(disease_id, "_consensus_reactome_enrichment.csv"))
      )

      plot_dt <- reactome_dt[
        order(p.adjust)
      ][
        1:min(15, .N)
      ]

      plot_dt[, pathway := stringr::str_wrap(Description, width = 45)]
      plot_dt[, neg_log10_padj := -log10(p.adjust)]
      plot_dt[, GeneRatio_num := {
        parts <- tstrsplit(GeneRatio, "/", fixed = TRUE)
        as.numeric(parts[[1]]) / as.numeric(parts[[2]])
      }]

      plot_dt[, pathway := forcats::fct_reorder(pathway, neg_log10_padj)]

      p_consensus_reactome <- ggplot(
        plot_dt,
        aes(
          x = neg_log10_padj,
          y = pathway,
          size = Count
        )
      ) +
        geom_point(
          colour = "#5B9BD5",
          alpha = 0.90
        ) +
        scale_size_continuous(
          name = "Gene count",
          range = c(3, 8)
        ) +
        labs(
          title = paste0("Consensus Reactome enrichment of ", disease_id, " feature-importance genes"),
          subtitle = "Top consensus genes from the final champion and comparator configurations",
          x = expression(-log[10]("adjusted p-value")),
          y = "Reactome pathway"
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5, colour = "black"),
          plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "black"),
          axis.text.y = element_text(size = 8.5, colour = "black"),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.title = element_text(face = "bold", colour = "black"),
          legend.title = element_text(face = "bold", colour = "black"),
          legend.text = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
          plot.margin = margin(15, 35, 15, 35)
        )

      print(p_consensus_reactome)

      ggsave(
        file.path(fig_dir, paste0(disease_id, "_consensus_reactome_enrichment.png")),
        plot = p_consensus_reactome,
        width = 260,
        height = 180,
        units = "mm",
        dpi = 600,
        bg = "white"
      )

      ggsave(
        file.path(fig_dir, paste0(disease_id, "_consensus_reactome_enrichment.pdf")),
        plot = p_consensus_reactome,
        width = 260,
        height = 180,
        units = "mm",
        device = cairo_pdf,
        bg = "white"
      )
    } else {
      warning("No Reactome pathways enriched for consensus genes.")
    }
  } else {
    warning("Too few mapped Entrez IDs for consensus Reactome enrichment.")
  }

  # ---------------------------------------------------------------------------
  # 9.2 Model-specific Reactome enrichment
  # ---------------------------------------------------------------------------

  top_n_genes <- 15

  model_gene_sets <- fi_summary[
    order(Configuration, -mean_loss)
  ][
    ,
    head(.SD, top_n_genes),
    by = Configuration
  ][
    ,
    .(
      genes = list(unique(variable))
    ),
    by = Configuration
  ]

  map_one_model <- function(config, genes) {
    mapped <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = genes,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    ) %>%
      as.data.table() %>%
      filter(!is.na(ENTREZID)) %>%
      distinct(SYMBOL, ENTREZID)

    mapped[, Configuration := config]
    mapped
  }

  mapped_all <- rbindlist(
    lapply(
      seq_len(nrow(model_gene_sets)),
      function(i) {
        map_one_model(
          model_gene_sets$Configuration[i],
          model_gene_sets$genes[[i]]
        )
      }
    ),
    fill = TRUE
  )

  fwrite(
    mapped_all,
    file.path(enrich_dir, paste0(disease_id, "_model_specific_top15_symbol_to_entrez.csv"))
  )

  reactome_by_model <- list()

  for (cfg in unique(mapped_all$Configuration)) {
    entrez_ids <- unique(mapped_all[Configuration == cfg, ENTREZID])

    if (length(entrez_ids) < 5) {
      warning("Skipping ", cfg, ": fewer than 5 mapped Entrez IDs.")
      next
    }

    enr <- ReactomePA::enrichPathway(
      gene = entrez_ids,
      organism = "human",
      pvalueCutoff = 0.10,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.20,
      readable = TRUE
    )

    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
      warning("No enriched Reactome pathways for: ", cfg)
      next
    }

    dt <- as.data.table(as.data.frame(enr))
    dt[, Configuration := cfg]

    reactome_by_model[[cfg]] <- dt
  }

  reactome_model_dt <- rbindlist(reactome_by_model, fill = TRUE)

  if (nrow(reactome_model_dt) > 0) {
    fwrite(
      reactome_model_dt,
      file.path(enrich_dir, paste0(disease_id, "_model_specific_top15_reactome_enrichment.csv"))
    )

    plot_model_dt <- reactome_model_dt[
      order(Configuration, p.adjust)
    ][
      ,
      head(.SD, 8),
      by = Configuration
    ]

    plot_model_dt[, pathway := stringr::str_wrap(Description, width = 42)]
    plot_model_dt[, neg_log10_padj := -log10(p.adjust)]
    plot_model_dt[, GeneRatio_num := {
      parts <- tstrsplit(GeneRatio, "/", fixed = TRUE)
      as.numeric(parts[[1]]) / as.numeric(parts[[2]])
    }]

    plot_model_dt[, Short_Config := as.character(Configuration)]

    if (requireNamespace("tidytext", quietly = TRUE)) {
      plot_model_dt[, pathway_ordered := tidytext::reorder_within(
        pathway,
        neg_log10_padj,
        Short_Config
      )]

      p_model_reactome <- ggplot(
        plot_model_dt,
        aes(
          x = neg_log10_padj,
          y = pathway_ordered,
          size = Count
        )
      ) +
        geom_point(
          colour = "#5B9BD5",
          alpha = 0.90
        ) +
        facet_wrap(
          ~ Short_Config,
          scales = "free_y",
          ncol = 2
        ) +
        tidytext::scale_y_reordered() +
        scale_size_continuous(
          name = "Gene count",
          range = c(2.5, 7.5)
        ) +
        labs(
          title = paste0("Model-specific Reactome enrichment of ", disease_id, " feature-importance genes"),
          subtitle = "Top Reactome pathways from the top 15 genes within each final configuration",
          x = expression(-log[10]("adjusted p-value")),
          y = "Reactome pathway"
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5, colour = "black"),
          plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "black"),
          axis.title = element_text(face = "bold", colour = "black"),
          axis.text.y = element_text(size = 7.5, colour = "black"),
          axis.text.x = element_text(size = 8, colour = "black"),
          strip.text = element_text(face = "bold", size = 8, colour = "black"),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          legend.title = element_text(face = "bold", colour = "black"),
          legend.text = element_text(colour = "black"),
          plot.margin = margin(10, 20, 10, 20)
        )

      print(p_model_reactome)

      ggsave(
        file.path(fig_dir, paste0(disease_id, "_model_specific_reactome_enrichment.png")),
        plot = p_model_reactome,
        width = 260,
        height = 210,
        units = "mm",
        dpi = 600,
        bg = "white"
      )

      ggsave(
        file.path(fig_dir, paste0(disease_id, "_model_specific_reactome_enrichment.pdf")),
        plot = p_model_reactome,
        width = 260,
        height = 210,
        units = "mm",
        device = cairo_pdf,
        bg = "white"
      )
    } else {
      warning("tidytext is not installed. Model-specific Reactome plot skipped.")
    }
  } else {
    warning("No model-specific Reactome enrichment results were obtained.")
  }
}

# ==============================================================================
# 10. Reproducibility log
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(summary_dir, paste0(disease_id, "_sessionInfo.txt"))
)

message("\n[done] Script 04 completed successfully.")
message("[done] Bundles: ", bundle_dir)
message("[done] Test outputs: ", test_dir)
message("[done] Feature importance: ", explain_dir)
message("[done] Reactome enrichment: ", enrich_dir)
message("[done] Figures: ", fig_dir)
