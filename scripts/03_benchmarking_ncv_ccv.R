# ==============================================================================
# TRACE-Omics
# Script 03: Benchmarking, nested cross-validation and cross-cohort validation
#
# Purpose:
#   This script evaluates transcriptomic classification pipelines across
#   harmonisation configurations, feature-selection frameworks, machine-learning
#   algorithms, class-balancing approaches, nested cross-validation and
#   cross-cohort validation.
#
# Required input structure:
#
#   data/processed/TB/tasks/feature_selection_tasks/
#   ├── train_TB_VST_only_mrmr.rds
#   ├── train_TB_VST_SVA_mrmr.rds
#   └── ...
#
#   data/processed/TB/tasks/balanced_tasks/
#   ├── train_task_TB_VST_SVA_mrmr_SMOTE.rds
#   ├── train_task_TB_VST_ComBatRef_HybridA_Upsample.rds
#   └── ...
#
# Expected mlr3 task requirement:
#   - Target column must be named: outcome
#   - Positive class must be: case
#   - Negative class must be: control
#   - For CCV, a non-feature column named cohort must be present
#
# Outputs:
#   - Baseline benchmark tables
#   - Baseline AUROC heatmap
#   - Prediction-level benchmark outputs
#   - ROC and calibration plots for selected configurations
#   - Nested CV fold-level and summary tables
#   - Hyperparameter sensitivity plots
#   - Outer CV predictions and metrics from best nested-CV parameters
#   - Permutation test table
#   - CCV prediction tables and cohort-level ROC/PR figures
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
  library(RColorBrewer)
  library(pROC)
  library(future)
  library(mlr3verse)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3tuning)
  library(mlr3mbo)
  library(mlr3extralearners)
  library(paradox)
  library(patchwork)
})

plan(sequential)
set.seed(7832)

lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

# ==============================================================================
# 1. User settings
# ==============================================================================

disease_id <- "TB"

task_dir <- file.path(
  "data", "processed", disease_id, "tasks", "feature_selection_tasks"
)

balanced_task_dir <- file.path(
  "data", "processed", disease_id, "tasks", "balanced_tasks"
)

base_out_dir <- file.path(
  "results", "03_benchmarking_ncv_ccv", disease_id
)

baseline_dir <- file.path(base_out_dir, "01_baseline_benchmark")
ncv_dir      <- file.path(base_out_dir, "02_nested_cv")
outercv_dir  <- file.path(base_out_dir, "03_outercv_from_bestparams")
ccv_dir      <- file.path(base_out_dir, "04_cross_cohort_validation")
fig_dir      <- file.path(base_out_dir, "figures")
summary_dir  <- file.path(base_out_dir, "summary")

dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ncv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outercv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ccv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

positive_class <- "case"
negative_class <- "control"

baseline_folds <- 5
baseline_repeats <- 2

outer_folds <- 5
inner_folds <- 5
n_evals <- 10

n_boot_metric <- 1000
n_boot_curve  <- 500

hc_order <- c(
  "ComBat-Seq + log1p",
  "VST only",
  "VST + SVA",
  "VST + ComBat",
  "VST + ComBat-Ref"
)

fsf_order <- c(
  "Boruta",
  "RF_Imp",
  "CMIM",
  "AUC",
  "mRMR",
  "Info_Gain",
  "ANOVA",
  "RF_Var",
  "Hybrid A",
  "Hybrid B"
)

mla_order <- c(
  "Featureless",
  "Naive Bayes",
  "Linear SVM",
  "RBF SVM",
  "C50",
  "Random Forest",
  "XGBoost",
  "LightGBM"
)

cba_order <- c("NBA", "Upsample", "SMOTE", "Downsample")

# ------------------------------------------------------------------------------
# Selected configurations for manuscript-style plots.
# Edit only these tables for TB/UC/RA figure-specific choices.
# ------------------------------------------------------------------------------

fig1_configs <- tribble(
  ~MLA,            ~HC,                  ~FSF,       ~label,
  "LightGBM",      "VST + ComBat-Ref",   "Hybrid A", "LightGBM",
  "XGBoost",       "VST + ComBat-Ref",   "Hybrid A", "XGBoost",
  "Random Forest", "VST + ComBat-Ref",   "Hybrid A", "Random Forest",
  "C50",           "VST + ComBat-Ref",   "Hybrid A", "C5.0",
  "RBF SVM",       "VST + ComBat-Ref",   "Hybrid A", "RBF SVM",
  "Linear SVM",    "VST + ComBat-Ref",   "Hybrid A", "Linear SVM",
  "Naive Bayes",   "VST + ComBat-Ref",   "Hybrid A", "Naive Bayes"
)

fig2_configs <- tribble(
  ~MLA,            ~HC,                    ~FSF,       ~label,
  "Random Forest", "VST + ComBat",         "Hybrid A", "VST + ComBat",
  "Random Forest", "VST + ComBat-Ref",     "Hybrid A", "VST + ComBat-Ref",
  "Random Forest", "VST + SVA",            "Hybrid A", "VST + SVA",
  "Random Forest", "VST only",             "Hybrid A", "VST only",
  "Random Forest", "ComBat-Seq + log1p",   "Hybrid A", "ComBat-Seq + log1p"
)

calibration_model_configs <- tribble(
  ~MLA,            ~HC,                  ~FSF,       ~label,
  "Random Forest", "VST + ComBat-Ref",   "Hybrid A", "Random Forest",
  "LightGBM",      "VST + ComBat-Ref",   "Hybrid A", "LightGBM",
  "XGBoost",       "VST + ComBat-Ref",   "Hybrid A", "XGBoost",
  "RBF SVM",       "VST + ComBat-Ref",   "Hybrid A", "RBF SVM",
  "Linear SVM",    "VST + ComBat-Ref",   "Hybrid A", "Linear SVM",
  "Naive Bayes",   "VST + ComBat-Ref",   "Hybrid A", "Naive Bayes",
  "C50",           "VST + ComBat-Ref",   "Hybrid A", "C50"
)

calibration_hc_configs <- fig2_configs

top4_configs_for_permutation <- c(
  "LightGBM | VST + ComBat-Ref | Info_Gain | Upsample",
  "Random Forest | VST + ComBat-Ref | Info_Gain | Upsample",
  "RBF SVM | VST + ComBat-Ref | Info_Gain | Upsample",
  "LightGBM | VST + ComBat | Hybrid A | Upsample"
)

ccv_winner <- list(
  MLA = "Random Forest",
  HC  = "VST + SVA",
  FSF = "mRMR",
  CBA = "SMOTE"
)

# ==============================================================================
# 2. General parsing and metric helpers
# ==============================================================================

clean_hc <- function(x) {
  case_when(
    str_detect(x, "VST_ComBatRef|ComBatRef")   ~ "VST + ComBat-Ref",
    str_detect(x, "VST_ComBat|ComBat")         ~ "VST + ComBat",
    str_detect(x, "VST_SVA|SVA")               ~ "VST + SVA",
    str_detect(x, "ComBatSeq_log1p")           ~ "ComBat-Seq + log1p",
    str_detect(x, "VST_only|VST only")         ~ "VST only",
    TRUE ~ as.character(x)
  )
}

clean_fsf <- function(x) {
  case_when(
    str_detect(x, "HybridA|Hybrid A")           ~ "Hybrid A",
    str_detect(x, "HybridB|Hybrid B")           ~ "Hybrid B",
    str_detect(x, "information_gain|Info_Gain") ~ "Info_Gain",
    str_detect(x, "boruta|Boruta")              ~ "Boruta",
    str_detect(x, "cmim|CMIM")                  ~ "CMIM",
    str_detect(x, "mrmr|mRMR")                  ~ "mRMR",
    str_detect(x, "anova|ANOVA")                ~ "ANOVA",
    str_detect(x, "auc|AUC")                    ~ "AUC",
    str_detect(x, "importance|RF_Imp")          ~ "RF_Imp",
    str_detect(x, "variance|RF_Var")            ~ "RF_Var",
    TRUE ~ NA_character_
  )
}

clean_cba <- function(x) {
  case_when(
    str_detect(x, "SMOTE$|SMOTE")           ~ "SMOTE",
    str_detect(x, "Upsample$|Upsample")     ~ "Upsample",
    str_detect(x, "Downsample$|Downsample") ~ "Downsample",
    str_detect(x, "NBA$|NBA")               ~ "NBA",
    TRUE ~ "NBA"
  )
}

clean_mla <- function(x) {
  case_when(
    str_detect(x, "Featureless")               ~ "Featureless",
    str_detect(x, "Naive")                     ~ "Naive Bayes",
    str_detect(x, "Linear_SVM|Linear SVM")     ~ "Linear SVM",
    str_detect(x, "RBF_SVM|RBF SVM")           ~ "RBF SVM",
    str_detect(x, "RandomForest|Random Forest") ~ "Random Forest",
    str_detect(x, "LightGBM")                  ~ "LightGBM",
    str_detect(x, "XGBoost")                   ~ "XGBoost",
    str_detect(x, "C50|C5")                    ~ "C50",
    TRUE ~ as.character(x)
  )
}

parse_task_id <- function(task_id) {
  tibble(
    task_id = task_id,
    HC  = clean_hc(task_id),
    FSF = clean_fsf(task_id),
    CBA = clean_cba(task_id)
  )
}

prob_col_name <- function(dt) {
  if ("prob.case" %in% names(dt)) return("prob.case")
  if ("prob_case" %in% names(dt)) return("prob_case")
  stop("Could not find probability column: expected prob.case or prob_case.", call. = FALSE)
}

calc_class_metrics <- function(truth, response, prob_case,
                               pos_class = positive_class,
                               neg_class = negative_class) {

  truth <- as.character(truth)
  response <- as.character(response)
  prob_case <- as.numeric(prob_case)

  tp <- sum(response == pos_class & truth == pos_class, na.rm = TRUE)
  tn <- sum(response == neg_class & truth == neg_class, na.rm = TRUE)
  fp <- sum(response == pos_class & truth == neg_class, na.rm = TRUE)
  fn <- sum(response == neg_class & truth == pos_class, na.rm = TRUE)

  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  bacc <- mean(c(sens, spec), na.rm = TRUE)

  y_true <- ifelse(truth == pos_class, 1, 0)
  brier <- mean((pmin(pmax(prob_case, 0), 1) - y_true)^2, na.rm = TRUE)

  c(
    BACC = bacc,
    Sens = sens,
    Spec = spec,
    Brier = brier
  )
}

bootstrap_metric_ci <- function(dt, truth_col, pred_col, prob_col,
                                B = n_boot_metric,
                                seed = 7832) {

  set.seed(seed)
  n <- nrow(dt)

  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    d <- dt[idx, ]

    calc_class_metrics(
      truth = d[[truth_col]],
      response = d[[pred_col]],
      prob_case = d[[prob_col]]
    )
  })

  data.table(
    BACC_l  = as.numeric(quantile(boot_mat["BACC", ], 0.025, na.rm = TRUE)),
    BACC_u  = as.numeric(quantile(boot_mat["BACC", ], 0.975, na.rm = TRUE)),
    Sens_l  = as.numeric(quantile(boot_mat["Sens", ], 0.025, na.rm = TRUE)),
    Sens_u  = as.numeric(quantile(boot_mat["Sens", ], 0.975, na.rm = TRUE)),
    Spec_l  = as.numeric(quantile(boot_mat["Spec", ], 0.025, na.rm = TRUE)),
    Spec_u  = as.numeric(quantile(boot_mat["Spec", ], 0.975, na.rm = TRUE)),
    Brier_l = as.numeric(quantile(boot_mat["Brier", ], 0.025, na.rm = TRUE)),
    Brier_u = as.numeric(quantile(boot_mat["Brier", ], 0.975, na.rm = TRUE))
  )
}

fmt_ci <- function(est, l, u) {
  sprintf("%.3f [%.3f--%.3f]", est, l, u)
}

bootstrap_mean_ci <- function(x, n_boot = 2000, conf = 0.95, seed = 7832) {
  set.seed(seed)

  x <- x[is.finite(x)]

  if (length(x) == 0) return(c(NA_real_, NA_real_))

  boot_means <- replicate(
    n_boot,
    mean(sample(x, size = length(x), replace = TRUE), na.rm = TRUE)
  )

  alpha <- (1 - conf) / 2
  quantile(boot_means, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
}

safe_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = path, compress = "gzip")
}

# ==============================================================================
# 3. Baseline benchmarking across HCs, FSFs and MLAs
# ==============================================================================

task_files <- list.files(
  task_dir,
  pattern = "^train_.*\\.rds$",
  full.names = TRUE
)

if (length(task_files) == 0) {
  stop("No baseline training task files found in: ", task_dir, call. = FALSE)
}

task_objects <- lapply(task_files, function(f) {
  tryCatch(readRDS(f), error = function(e) NULL)
})

is_task <- vapply(task_objects, inherits, logical(1), what = "TaskClassif")

task_files <- task_files[is_task]
tasks <- task_objects[is_task]
names(tasks) <- basename(task_files)

message("[baseline] valid tasks: ", length(tasks))

learners <- list(
  lrn("classif.featureless", predict_type = "prob", id = "Featureless"),
  lrn("classif.C50", predict_type = "prob", id = "C50"),
  lrn("classif.ranger", predict_type = "prob", id = "Random Forest"),
  lrn("classif.xgboost", predict_type = "prob", id = "XGBoost"),
  lrn("classif.lightgbm", predict_type = "prob", id = "LightGBM"),
  lrn("classif.naive_bayes", predict_type = "prob", id = "Naive Bayes"),
  lrn("classif.svm", predict_type = "prob", kernel = "linear", id = "Linear SVM"),
  lrn("classif.svm", predict_type = "prob", kernel = "radial", id = "RBF SVM")
)

baseline_resampling <- rsmp(
  "repeated_cv",
  folds = baseline_folds,
  repeats = baseline_repeats
)

baseline_design <- benchmark_grid(
  tasks = tasks,
  learners = learners,
  resamplings = baseline_resampling
)

bmr <- benchmark(baseline_design)

score_dt <- as.data.table(
  bmr$score(msrs(c("classif.auc", "classif.bacc", "classif.bbrier")))
)

agg_dt <- as.data.table(
  bmr$aggregate(msrs(c("classif.auc", "classif.bacc", "classif.bbrier")))
)

agg_dt_clean <- copy(agg_dt)
if ("resample_result" %in% names(agg_dt_clean)) {
  agg_dt_clean[, resample_result := NULL]
}

fwrite(score_dt, file.path(baseline_dir, "baseline_fold_level_scores.csv"))
fwrite(agg_dt_clean, file.path(baseline_dir, "baseline_aggregated_scores.csv"))

saveRDS(bmr, file.path(baseline_dir, "baseline_benchmark_result.rds"))
saveRDS(agg_dt, file.path(baseline_dir, "baseline_agg_dt_full.rds"))
saveRDS(agg_dt_clean, file.path(baseline_dir, "baseline_agg_dt_clean.rds"))

# ==============================================================================
# 4. Derive benchmark labels and baseline heatmap
# ==============================================================================

agg_dt_labelled <- agg_dt_clean |>
  as_tibble() |>
  rename(MLA = learner_id) |>
  mutate(
    MLA = clean_mla(MLA),
    HC = clean_hc(task_id),
    FSF = clean_fsf(task_id),
    Mean_AUROC = classif.auc
  ) |>
  filter(!is.na(HC), !is.na(FSF), !is.na(MLA), !is.na(Mean_AUROC)) |>
  mutate(
    HC = factor(HC, levels = hc_order),
    FSF = factor(FSF, levels = fsf_order),
    MLA = factor(MLA, levels = mla_order)
  )

fwrite(
  as.data.table(agg_dt_labelled),
  file.path(baseline_dir, "baseline_aggregated_scores_labelled.csv")
)

heat_main <- agg_dt_labelled |>
  filter(MLA %in% mla_order)

p_heat <- ggplot(heat_main, aes(x = FSF, y = HC, fill = Mean_AUROC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", Mean_AUROC)), size = 2.3, colour = "black") +
  facet_wrap(~ MLA, ncol = 2) +
  scale_fill_viridis_c(
    name = "AUROC",
    limits = c(0.5, 1),
    oob = scales::squish
  ) +
  labs(
    title = paste0(disease_id, ": baseline benchmarking across HC-FSF-MLA configurations"),
    x = "Feature-selection framework",
    y = "Harmonisation configuration"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, paste0(disease_id, "_baseline_AUROC_heatmap.png")),
  p_heat,
  width = 240,
  height = 200,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_baseline_AUROC_heatmap.pdf")),
  p_heat,
  width = 240,
  height = 200,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 5. Extract baseline prediction-level outputs
# ==============================================================================

extract_prediction_dt <- function(bmr) {
  rows <- list()
  rr_list <- bmr$resample_results$resample_result

  for (i in seq_along(rr_list)) {
    rr <- rr_list[[i]]

    task_id <- bmr$resample_results$task_id[i]
    learner_id <- bmr$resample_results$learner_id[i]
    resampling_id <- bmr$resample_results$resampling_id[i]

    preds <- rr$predictions()

    pred_rows <- imap(preds, function(pred, iteration_id) {
      dt <- as.data.table(pred)

      dt[, `:=`(
        task_id = task_id,
        learner_id = learner_id,
        resampling_id = resampling_id,
        iteration = iteration_id,
        MLA = clean_mla(learner_id),
        HC = clean_hc(task_id),
        FSF = clean_fsf(task_id),
        CBA = clean_cba(task_id)
      )]

      dt
    })

    rows[[i]] <- rbindlist(pred_rows, fill = TRUE)
  }

  rbindlist(rows, fill = TRUE)
}

pred_dt <- extract_prediction_dt(bmr)

pcol <- prob_col_name(pred_dt)
if (pcol == "prob.case") {
  pred_dt[, prob_case := `prob.case`]
}

pred_dt[, Configuration := paste(MLA, HC, FSF, CBA, sep = " | ")]

fwrite(
  pred_dt,
  file.path(baseline_dir, "baseline_prediction_level.csv")
)

# ==============================================================================
# 6. Baseline scalar metrics with CIs
# ==============================================================================

calc_metrics_by_config <- function(dt) {
  truth <- as.character(dt$truth)
  score <- as.numeric(dt$prob_case)
  pred  <- as.character(dt$response)

  if (length(unique(truth)) < 2) {
    return(data.table(
      AUROC = NA_real_, AUROC_l = NA_real_, AUROC_u = NA_real_,
      BACC = NA_real_, BACC_l = NA_real_, BACC_u = NA_real_,
      Sens = NA_real_, Sens_l = NA_real_, Sens_u = NA_real_,
      Spec = NA_real_, Spec_l = NA_real_, Spec_u = NA_real_,
      Brier = NA_real_, Brier_l = NA_real_, Brier_u = NA_real_
    ))
  }

  roc_obj <- pROC::roc(
    response = truth,
    predictor = score,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))

  point_metrics <- calc_class_metrics(
    truth = truth,
    response = pred,
    prob_case = score
  )

  boot_ci <- bootstrap_metric_ci(
    dt = as.data.table(dt),
    truth_col = "truth",
    pred_col = "response",
    prob_col = "prob_case",
    B = n_boot_metric
  )

  data.table(
    AUROC = as.numeric(pROC::auc(roc_obj)),
    AUROC_l = auc_ci[1],
    AUROC_u = auc_ci[3],
    BACC = as.numeric(point_metrics["BACC"]),
    BACC_l = boot_ci$BACC_l,
    BACC_u = boot_ci$BACC_u,
    Sens = as.numeric(point_metrics["Sens"]),
    Sens_l = boot_ci$Sens_l,
    Sens_u = boot_ci$Sens_u,
    Spec = as.numeric(point_metrics["Spec"]),
    Spec_l = boot_ci$Spec_l,
    Spec_u = boot_ci$Spec_u,
    Brier = as.numeric(point_metrics["Brier"]),
    Brier_l = boot_ci$Brier_l,
    Brier_u = boot_ci$Brier_u
  )
}

baseline_metrics_ci <- pred_dt[
  ,
  calc_metrics_by_config(.SD),
  by = .(MLA, HC, FSF, CBA)
]

baseline_metrics_pretty <- baseline_metrics_ci |>
  as_tibble() |>
  mutate(
    AUROC_CI = fmt_ci(AUROC, AUROC_l, AUROC_u),
    BACC_CI  = fmt_ci(BACC, BACC_l, BACC_u),
    Sens_CI  = fmt_ci(Sens, Sens_l, Sens_u),
    Spec_CI  = fmt_ci(Spec, Spec_l, Spec_u),
    Brier_CI = fmt_ci(Brier, Brier_l, Brier_u),
    Configuration = paste(MLA, HC, FSF, CBA, sep = " | ")
  ) |>
  arrange(desc(AUROC), desc(BACC))

fwrite(
  as.data.table(baseline_metrics_ci),
  file.path(baseline_dir, "baseline_metrics_with_ci.csv")
)

fwrite(
  as.data.table(baseline_metrics_pretty),
  file.path(baseline_dir, "baseline_metrics_with_ci_pretty.csv")
)

# ==============================================================================
# 7. ROC helper and selected ROC plots
# ==============================================================================

plot_selected_roc <- function(pred_dt, config_tbl, title, subtitle,
                              boot_n = n_boot_curve,
                              ci_grid = seq(0, 1, length.out = 100)) {

  plot_dt <- as_tibble(pred_dt) |>
    inner_join(config_tbl, by = c("MLA", "HC", "FSF")) |>
    mutate(
      prob_case = as.numeric(prob_case),
      label = as.character(label)
    ) |>
    filter(!is.na(prob_case)) |>
    group_by(label) |>
    filter(n_distinct(truth) == 2) |>
    ungroup()

  if (nrow(plot_dt) == 0) {
    stop("No predictions matched selected ROC configurations.", call. = FALSE)
  }

  roc_results <- plot_dt |>
    group_by(label) |>
    group_split() |>
    map(function(df) {
      lab <- unique(df$label)

      roc_obj <- pROC::roc(
        df$truth,
        df$prob_case,
        levels = c(negative_class, positive_class),
        direction = "<",
        quiet = TRUE
      )

      roc_coords <- pROC::coords(
        roc_obj,
        x = ci_grid,
        input = "specificity",
        ret = c("specificity", "sensitivity"),
        transpose = FALSE
      ) |>
        as.data.frame() |>
        mutate(
          fpr = 1 - specificity,
          label = lab
        )

      auc_ci <- pROC::ci.auc(roc_obj, method = "delong")

      ci_se <- pROC::ci.se(
        roc_obj,
        specificities = ci_grid,
        boot.n = boot_n
      )

      ci_df <- as.data.frame(ci_se)
      ci_df$specificity <- as.numeric(rownames(ci_df))

      ci_df <- ci_df |>
        setNames(c("sens_lower", "sensitivity", "sens_upper", "specificity")) |>
        mutate(
          fpr = 1 - specificity,
          label = lab
        )

      list(
        roc_df = roc_coords,
        ci_df = ci_df,
        auc_ci = auc_ci
      )
    })

  roc_df_all <- map_dfr(roc_results, "roc_df")
  ci_df_all  <- map_dfr(roc_results, "ci_df")

  legend_tbl <- map_dfr(roc_results, function(x) {
    tibble(
      label = unique(x$roc_df$label),
      AUROC = as.numeric(x$auc_ci[2]),
      low = as.numeric(x$auc_ci[1]),
      high = as.numeric(x$auc_ci[3])
    )
  }) |>
    arrange(desc(AUROC)) |>
    mutate(
      legend_label = sprintf(
        "%s (AUROC %.3f [%.3f--%.3f])",
        label, AUROC, low, high
      )
    )

  roc_df_all <- roc_df_all |>
    left_join(legend_tbl, by = "label") |>
    mutate(legend_label = factor(legend_label, levels = legend_tbl$legend_label))

  ci_df_all <- ci_df_all |>
    left_join(legend_tbl, by = "label") |>
    mutate(legend_label = factor(legend_label, levels = legend_tbl$legend_label))

  pal <- setNames(
    RColorBrewer::brewer.pal(max(3, nrow(legend_tbl)), "Set1")[seq_len(nrow(legend_tbl))],
    legend_tbl$legend_label
  )

  ggplot() +
    geom_ribbon(
      data = ci_df_all,
      aes(x = fpr, ymin = sens_lower, ymax = sens_upper, fill = legend_label),
      alpha = 0.12,
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
      aes(x = fpr, y = sensitivity, colour = legend_label),
      linewidth = 1.15,
      lineend = "round"
    ) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    coord_equal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      axis.text = element_text(colour = "black", size = 9),
      axis.title = element_text(colour = "black", face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11.5, hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 8.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
    )
}

fig1_roc <- plot_selected_roc(
  pred_dt,
  fig1_configs,
  title = "Performance Ceiling",
  subtitle = "AUROC for selected MLA-HC-FSF configurations"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_figure1_performance_ceiling_ROC.png")),
  fig1_roc,
  width = 8.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_figure1_performance_ceiling_ROC.pdf")),
  fig1_roc,
  width = 8.5,
  height = 6.5,
  device = cairo_pdf,
  bg = "white"
)

fig2_roc <- plot_selected_roc(
  pred_dt,
  fig2_configs,
  title = "Effect of Harmonisation Configurations",
  subtitle = "Random Forest model with Hybrid A feature set"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_figure2_hc_effect_ROC.png")),
  fig2_roc,
  width = 8.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_figure2_hc_effect_ROC.pdf")),
  fig2_roc,
  width = 8.5,
  height = 6.5,
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 8. Calibration helper and selected calibration plots
# ==============================================================================

compute_calibration_curve <- function(y, p, n_bins = 10) {
  tibble(y = y, p = pmin(pmax(p, 0), 1)) |>
    mutate(bin = ntile(p, n_bins)) |>
    group_by(bin) |>
    summarise(
      mean_p = mean(p, na.rm = TRUE),
      mean_y = mean(y, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) |>
    arrange(mean_p)
}

calc_brier <- function(y, p) {
  mean((pmin(pmax(p, 0), 1) - y)^2, na.rm = TRUE)
}

calc_ece <- function(cal_df) {
  sum((cal_df$n / sum(cal_df$n)) * abs(cal_df$mean_y - cal_df$mean_p))
}

bootstrap_calibration_band <- function(y, p,
                                       grid = seq(0, 1, by = 0.01),
                                       B = n_boot_curve,
                                       n_bins = 10) {

  n <- length(y)

  if (n < 10 || length(unique(y)) < 2) {
    return(
      tibble(
        mean_p = grid,
        obs_low = NA_real_,
        obs_mid = NA_real_,
        obs_high = NA_real_
      )
    )
  }

  mats <- replicate(B, {
    idx <- sample.int(n, size = n, replace = TRUE)
    yy <- y[idx]
    pp <- p[idx]

    if (length(unique(yy)) < 2) {
      return(rep(NA_real_, length(grid)))
    }

    cal <- compute_calibration_curve(yy, pp, n_bins = n_bins)

    approx(
      x = cal$mean_p,
      y = cal$mean_y,
      xout = grid,
      method = "linear",
      rule = 2,
      ties = "ordered"
    )$y
  })

  if (is.vector(mats)) {
    mats <- matrix(mats, nrow = length(grid))
  }

  tibble(
    mean_p = grid,
    obs_low = apply(mats, 1, quantile, probs = 0.025, na.rm = TRUE),
    obs_mid = apply(mats, 1, median, na.rm = TRUE),
    obs_high = apply(mats, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

plot_selected_calibration <- function(pred_dt, config_tbl, title, subtitle,
                                      n_bins = 10,
                                      boot_n = n_boot_curve,
                                      cal_grid = seq(0, 1, by = 0.01)) {

  plot_dt <- as_tibble(pred_dt) |>
    inner_join(config_tbl, by = c("MLA", "HC", "FSF")) |>
    mutate(
      y = ifelse(as.character(truth) == positive_class, 1, 0),
      p = as.numeric(prob_case),
      label = as.character(label)
    ) |>
    filter(!is.na(p)) |>
    group_by(label) |>
    filter(n_distinct(y) == 2) |>
    ungroup()

  if (nrow(plot_dt) == 0) {
    stop("No predictions matched selected calibration configurations.", call. = FALSE)
  }

  cal_list <- plot_dt |>
    group_by(label) |>
    group_split() |>
    map(function(d) {

      cal_df <- compute_calibration_curve(
        y = d$y,
        p = d$p,
        n_bins = n_bins
      ) |>
        mutate(label = unique(d$label))

      rib_df <- bootstrap_calibration_band(
        y = d$y,
        p = d$p,
        grid = cal_grid,
        B = boot_n,
        n_bins = n_bins
      ) |>
        mutate(label = unique(d$label))

      brier <- calc_brier(d$y, d$p)
      ece <- calc_ece(cal_df)

      tibble(
        label = unique(d$label),
        Brier = brier,
        ECE = ece,
        cal_df = list(cal_df),
        rib_df = list(rib_df)
      )
    })

  cal_df <- bind_rows(map(cal_list, "cal_df"))
  ribbon_df <- bind_rows(map(cal_list, "rib_df"))

  legend_tbl <- bind_rows(map(cal_list, ~ select(.x, label, Brier, ECE))) |>
    distinct() |>
    arrange(Brier) |>
    mutate(
      legend_label = sprintf(
        "%s (Brier %.3f; ECE %.3f)",
        label, Brier, ECE
      )
    )

  cal_df <- cal_df |>
    left_join(legend_tbl, by = "label") |>
    mutate(legend_label = factor(legend_label, levels = legend_tbl$legend_label))

  ribbon_df <- ribbon_df |>
    left_join(legend_tbl, by = "label") |>
    mutate(legend_label = factor(legend_label, levels = legend_tbl$legend_label))

  pal <- setNames(
    RColorBrewer::brewer.pal(max(3, nrow(legend_tbl)), "Set1")[seq_len(nrow(legend_tbl))],
    legend_tbl$legend_label
  )

  ggplot() +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.7
    ) +
    geom_ribbon(
      data = ribbon_df,
      aes(x = mean_p, ymin = obs_low, ymax = obs_high, fill = legend_label),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(
      data = cal_df,
      aes(x = mean_p, y = mean_y, colour = legend_label),
      linewidth = 1.25,
      lineend = "round"
    ) +
    geom_point(
      data = cal_df,
      aes(x = mean_p, y = mean_y, colour = legend_label),
      size = 2.3,
      alpha = 0.95
    ) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    coord_fixed() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Predicted probability (case)",
      y = "Observed fraction of cases"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, colour = "black"),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      legend.text = element_text(size = 9.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

cal_model <- plot_selected_calibration(
  pred_dt = pred_dt,
  config_tbl = calibration_model_configs,
  title = "Calibration comparison: model capacity under optimal representation",
  subtitle = "Selected HC-FSF representation used to isolate MLA capacity"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_model_capacity_calibration.png")),
  cal_model,
  width = 8.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_model_capacity_calibration.pdf")),
  cal_model,
  width = 8.5,
  height = 6.5,
  bg = "white"
)

cal_hc <- plot_selected_calibration(
  pred_dt = pred_dt,
  config_tbl = calibration_hc_configs,
  title = "Effect of Harmonisation Configurations",
  subtitle = "Random Forest model with Hybrid A feature set"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_hc_effect_calibration.png")),
  cal_hc,
  width = 8.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, paste0(disease_id, "_hc_effect_calibration.pdf")),
  cal_hc,
  width = 8.5,
  height = 6.5,
  bg = "white"
)

# ==============================================================================
# 9. Nested CV tuning
# ==============================================================================

learners_param_space <- list(
  RandomForest = list(
    label = "Random Forest",
    learner_id = "classif.ranger",
    param_set = ps(
      classif.ranger.num.trees     = p_int(200, 300),
      classif.ranger.mtry.ratio    = p_dbl(0.10, 0.5),
      classif.ranger.min.node.size = p_int(1, 8)
    ),
    learner_params = list(
      num.threads = 1,
      importance = "none",
      save.memory = TRUE
    )
  ),

  LightGBM = list(
    label = "LightGBM",
    learner_id = "classif.lightgbm",
    param_set = ps(
      classif.lightgbm.learning_rate    = p_dbl(0.01, 0.25),
      classif.lightgbm.num_leaves       = p_int(20, 100),
      classif.lightgbm.feature_fraction = p_dbl(0.6, 1.0)
    ),
    learner_params = list()
  ),

  RBF_SVM = list(
    label = "RBF SVM",
    learner_id = "classif.svm",
    param_set = ps(
      classif.svm.cost  = p_dbl(1e-2, 1e2, logscale = TRUE),
      classif.svm.gamma = p_dbl(1e-5, 2, logscale = TRUE)
    ),
    learner_params = list(
      type = "C-classification",
      kernel = "radial"
    )
  )
)

outer_cv <- rsmp("cv", folds = outer_folds)
inner_cv <- rsmp("cv", folds = inner_folds)

tuner_mbo <- tnr("mbo")
terminator <- trm("evals", n_evals = n_evals)
measure <- msr("classif.bacc")

build_autotuner <- function(learner_id, param_set, learner_params = list()) {

  base_args <- c(
    list(
      learner_id,
      predict_type = "prob"
    ),
    learner_params
  )

  learner_graph <- as_learner(
    po("removeconstants") %>>%
      do.call(lrn, base_args)
  )

  AutoTuner$new(
    learner = learner_graph,
    resampling = inner_cv$clone(),
    measure = measure,
    search_space = param_set,
    tuner = tuner_mbo,
    terminator = terminator,
    store_models = TRUE
  )
}

ncv_task_files <- list.files(
  balanced_task_dir,
  pattern = "^train_task_.*\\.rds$",
  full.names = TRUE
)

if (length(ncv_task_files) == 0) {
  warning("No balanced tasks found for nested CV in: ", balanced_task_dir)
}

names(ncv_task_files) <- basename(ncv_task_files) |>
  str_remove("\\.rds$") |>
  str_remove("^train_task_")

fold_rows <- list()
summary_rows <- list()
runtime_param_rows <- list()

for (task_name in names(ncv_task_files)) {

  task_raw <- readRDS(ncv_task_files[[task_name]])
  task_dt <- as.data.table(task_raw$data(cols = task_raw$col_names))

  task <- TaskClassif$new(
    id = task_name,
    backend = task_dt,
    target = "outcome",
    positive = positive_class
  )

  task$col_roles$stratum <- "outcome"
  task$col_roles$group <- character(0)

  meta <- parse_task_id(task_name)

  if (length(task$feature_names) < 2) {
    warning("Skipping nested CV task with fewer than 2 features: ", task_name)
    next
  }

  for (learner_name in names(learners_param_space)) {

    spec <- learners_param_space[[learner_name]]

    message("\n[nCV] Task: ", task_name, " | Learner: ", spec$label)

    at <- build_autotuner(
      learner_id = spec$learner_id,
      param_set = spec$param_set,
      learner_params = spec$learner_params
    )

    rr <- resample(
      task = task,
      learner = at,
      resampling = outer_cv$clone(),
      store_models = TRUE
    )

    rr_path <- file.path(
      ncv_dir,
      paste0("nestedcv_", learner_name, "_task_", task_name, ".rds")
    )

    safe_save_rds(rr, rr_path)

    fold_scores <- rr$score(measures = list(msr("classif.bacc")))
    bacc_col <- grep("classif\\.bacc$", names(fold_scores), value = TRUE)

    fold_dt <- as.data.table(fold_scores)
    fold_dt[, learner := spec$label]
    fold_dt[, learner_key := learner_name]
    fold_dt[, task_name := task_name]
    fold_dt[, HC := meta$HC]
    fold_dt[, FSF := meta$FSF]
    fold_dt[, CBA := meta$CBA]
    fold_dt[, fold_bacc := get(bacc_col)]

    fold_rows[[paste(task_name, learner_name, sep = "_")]] <- fold_dt[, .(
      task_name, HC, FSF, CBA, learner, learner_key, iteration, fold_bacc
    )]

    x <- fold_dt$fold_bacc
    ci <- bootstrap_mean_ci(x, n_boot = 2000, conf = 0.95, seed = 7832)

    summary_rows[[paste(task_name, learner_name, sep = "_")]] <- tibble(
      Model = spec$label,
      Learner_key = learner_name,
      HC = meta$HC,
      CBA = meta$CBA,
      FSF = meta$FSF,
      Mean_BACC = mean(x, na.rm = TRUE),
      Median_BACC = median(x, na.rm = TRUE),
      SD_BACC = sd(x, na.rm = TRUE),
      CI_low = ci[1],
      CI_high = ci[2],
      `95%_CI` = sprintf("(%.3f-%.3f)", ci[1], ci[2]),
      N_folds = length(x),
      RDS_path = rr_path
    )

    for (i in seq_along(rr$learners)) {
      learner_i <- rr$learners[[i]]
      archive_best <- tryCatch(learner_i$archive$best(), error = function(e) NULL)

      if (!is.null(archive_best) && length(archive_best$x_domain) > 0) {
        params <- archive_best$x_domain[[1]]

        runtime_param_rows[[length(runtime_param_rows) + 1]] <- tibble(
          task_name = task_name,
          HC = meta$HC,
          FSF = meta$FSF,
          CBA = meta$CBA,
          MLA = spec$label,
          fold = i,
          BACC = x[i],
          learning_rate = as.numeric(params[["classif.lightgbm.learning_rate"]] %||% NA),
          num_leaves = as.numeric(params[["classif.lightgbm.num_leaves"]] %||% NA),
          feature_fraction = as.numeric(params[["classif.lightgbm.feature_fraction"]] %||% NA),
          mtry.ratio = as.numeric(params[["classif.ranger.mtry.ratio"]] %||% NA),
          num.trees = as.numeric(params[["classif.ranger.num.trees"]] %||% NA),
          min.node.size = as.numeric(params[["classif.ranger.min.node.size"]] %||% NA),
          cost = as.numeric(params[["classif.svm.cost"]] %||% NA),
          gamma = as.numeric(params[["classif.svm.gamma"]] %||% NA)
        )
      }
    }
  }
}

fold_tbl <- bind_rows(fold_rows)

summary_tbl <- bind_rows(summary_rows) |>
  arrange(desc(Mean_BACC), desc(Median_BACC))

runtime_param_tbl <- bind_rows(runtime_param_rows)

write_csv(
  fold_tbl,
  file.path(summary_dir, "nestedcv_fold_level_results.csv")
)

write_csv(
  summary_tbl,
  file.path(summary_dir, "nestedcv_summary_results.csv")
)

write_csv(
  runtime_param_tbl,
  file.path(summary_dir, "nestedcv_runtime_hyperparameters_fold_level.csv")
)

# ==============================================================================
# 10. Hyperparameter sensitivity
# ==============================================================================

plot_hp_trend_pub <- function(data, xvar, xlab, model_label) {

  data_plot <- data |>
    filter(!is.na(.data[[xvar]]), !is.na(BACC), !is.na(CBA)) |>
    mutate(CBA = factor(CBA, levels = cba_order))

  ggplot(data_plot, aes(x = .data[[xvar]], y = BACC)) +
    geom_point(colour = "grey35", alpha = 0.50, size = 2.0) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 4),
      se = TRUE,
      colour = "black",
      fill = "grey72",
      linewidth = 1.4
    ) +
    facet_wrap(~ CBA, nrow = 1, drop = FALSE) +
    scale_y_continuous(
      limits = c(0.40, 1.00),
      breaks = seq(0.40, 1.00, by = 0.10),
      labels = number_format(accuracy = 0.01)
    ) +
    labs(
      title = model_label,
      x = xlab,
      y = "BACC"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 17, face = "bold", colour = "black", hjust = 0.5),
      axis.title = element_text(size = 12.5, face = "bold", colour = "black"),
      axis.text = element_text(size = 10.5, colour = "black"),
      strip.text = element_text(size = 12.5, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

if (nrow(runtime_param_tbl) > 0) {
  runtime_param_tbl_clean <- runtime_param_tbl |>
    mutate(
      CBA = factor(CBA, levels = cba_order),
      MLA = factor(MLA, levels = c("RBF SVM", "Random Forest", "LightGBM")),
      BACC = as.numeric(BACC)
    )

  hp_plots <- list()

  if (any(runtime_param_tbl_clean$MLA == "LightGBM", na.rm = TRUE)) {
    hp_plots$lgb <- runtime_param_tbl_clean |>
      filter(MLA == "LightGBM") |>
      plot_hp_trend_pub("learning_rate", "Learning rate", "LightGBM")
  }

  if (any(runtime_param_tbl_clean$MLA == "Random Forest", na.rm = TRUE)) {
    hp_plots$rf <- runtime_param_tbl_clean |>
      filter(MLA == "Random Forest") |>
      plot_hp_trend_pub("mtry.ratio", "mtry ratio", "Random Forest")
  }

  if (any(runtime_param_tbl_clean$MLA == "RBF SVM", na.rm = TRUE)) {
    hp_plots$svm <- runtime_param_tbl_clean |>
      filter(MLA == "RBF SVM") |>
      plot_hp_trend_pub("cost", "Cost (C)", "RBF SVM")
  }

  if (length(hp_plots) > 0) {
    hp_composite <- wrap_plots(hp_plots, ncol = 1) +
      plot_annotation(
        title = paste0("Hyperparameter sensitivity across selected MLAs for ", disease_id, " classification"),
        subtitle = "Fold-level BACC trends across class-balancing approaches during nested cross-validation",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 13, hjust = 0.5)
        )
      )

    ggsave(
      filename = file.path(fig_dir, paste0(disease_id, "_hyperparameter_sensitivity_selected_MLAs.png")),
      plot = hp_composite,
      width = 15,
      height = 9.5,
      dpi = 600,
      bg = "white"
    )

    ggsave(
      filename = file.path(fig_dir, paste0(disease_id, "_hyperparameter_sensitivity_selected_MLAs.pdf")),
      plot = hp_composite,
      device = cairo_pdf,
      width = 15,
      height = 9.5,
      units = "in",
      bg = "white"
    )
  }
}

# ==============================================================================
# 11. Outer CV from best nested-CV parameters
# ==============================================================================

extract_best_params <- function(rr, score_col = "classif.bacc") {

  learner_obj <- rr$learners[[1]]

  archive_data <- tryCatch(
    learner_obj$archive$data,
    error = function(e) NULL
  )

  if (is.null(archive_data)) {
    stop("No tuning archive found.", call. = FALSE)
  }

  if (!score_col %in% names(archive_data)) {
    stop("Could not find ", score_col, " in archive$data.", call. = FALSE)
  }

  best_idx <- which.max(archive_data[[score_col]])

  if ("x_domain" %in% names(archive_data)) {
    return(archive_data$x_domain[[best_idx]])
  }

  as.list(archive_data[best_idx, ])
}

build_fixed_learner <- function(model, best_params) {

  base_graph <- po("encode", method = "one-hot") %>>%
    po("removeconstants")

  if (model == "LightGBM") {
    return(as_learner(
      base_graph %>>%
        lrn(
          "classif.lightgbm",
          predict_type = "prob",
          learning_rate = as.numeric(best_params[["classif.lightgbm.learning_rate"]]),
          num_leaves = as.integer(best_params[["classif.lightgbm.num_leaves"]]),
          feature_fraction = as.numeric(best_params[["classif.lightgbm.feature_fraction"]])
        )
    ))
  }

  if (model == "Random Forest") {
    return(as_learner(
      base_graph %>>%
        lrn(
          "classif.ranger",
          predict_type = "prob",
          num.threads = 1,
          importance = "none",
          save.memory = TRUE,
          num.trees = as.integer(best_params[["classif.ranger.num.trees"]]),
          mtry.ratio = as.numeric(best_params[["classif.ranger.mtry.ratio"]]),
          min.node.size = as.integer(best_params[["classif.ranger.min.node.size"]])
        )
    ))
  }

  if (model == "RBF SVM") {
    return(as_learner(
      base_graph %>>%
        lrn(
          "classif.svm",
          predict_type = "prob",
          type = "C-classification",
          kernel = "radial",
          cost = as.numeric(best_params[["classif.svm.cost"]]),
          gamma = as.numeric(best_params[["classif.svm.gamma"]])
        )
    ))
  }

  stop("Unsupported model: ", model, call. = FALSE)
}

nestedcv_files <- list.files(
  ncv_dir,
  pattern = "^nestedcv_.*\\.rds$",
  full.names = TRUE
)

outercv_rows <- list()

for (f in nestedcv_files) {

  rr <- readRDS(f)
  fname <- basename(f)

  clean_name <- fname |>
    str_remove("^nestedcv_") |>
    str_remove("\\.rds$")

  learner_key <- str_extract(clean_name, "^(RandomForest|LightGBM|RBF_SVM)")
  task_name <- clean_name |>
    str_remove("^(RandomForest|LightGBM|RBF_SVM)_task_")

  model_label <- learners_param_space[[learner_key]]$label
  task_path <- file.path(balanced_task_dir, paste0("train_task_", task_name, ".rds"))

  if (!file.exists(task_path)) {
    warning("Could not find task for outer CV: ", task_path)
    next
  }

  task_raw <- readRDS(task_path)
  task_dt <- as.data.table(task_raw$data(cols = task_raw$col_names))

  task <- TaskClassif$new(
    id = task_name,
    backend = task_dt,
    target = "outcome",
    positive = positive_class
  )

  task$col_roles$stratum <- "outcome"
  task$col_roles$group <- character(0)

  best_params <- extract_best_params(rr)
  fixed_learner <- build_fixed_learner(model_label, best_params)

  outer_rr <- resample(
    task = task,
    learner = fixed_learner,
    resampling = outer_cv$clone(),
    store_models = TRUE
  )

  outer_path <- file.path(
    outercv_dir,
    paste0("outercv_nestedcv_", learner_key, "_task_", task_name, ".rds")
  )

  safe_save_rds(outer_rr, outer_path)

  outercv_rows[[length(outercv_rows) + 1]] <- tibble(
    File = basename(outer_path),
    Learner_key = learner_key,
    MLA = model_label,
    task_name = task_name,
    RDS_path = outer_path
  )
}

outercv_index <- bind_rows(outercv_rows)

write_csv(
  outercv_index,
  file.path(outercv_dir, "outercv_result_index.csv")
)

# ==============================================================================
# 12. Extract outer CV predictions and metrics
# ==============================================================================

outer_rds_files <- list.files(
  outercv_dir,
  pattern = "^outercv_nestedcv_.*\\.rds$",
  full.names = TRUE
)

prediction_rows <- list()

for (f in outer_rds_files) {

  obj <- readRDS(f)

  if (!inherits(obj, "ResampleResult")) next

  pred <- obj$prediction()
  pred_dt_tmp <- as.data.table(pred)

  fname <- basename(f)

  config <- fname |>
    str_remove("^outercv_nestedcv_") |>
    str_remove("\\.rds$")

  MLA <- clean_mla(config)
  HC <- clean_hc(config)
  FSF <- clean_fsf(config)
  CBA <- clean_cba(config)

  pred_dt_tmp[, `:=`(
    MLA = MLA,
    HC = HC,
    FSF = FSF,
    CBA = CBA,
    Configuration = paste(MLA, HC, FSF, CBA, sep = " | "),
    source_file = fname
  )]

  prediction_rows[[length(prediction_rows) + 1]] <- pred_dt_tmp
}

outer_pred_dt <- rbindlist(prediction_rows, fill = TRUE)

if (nrow(outer_pred_dt) > 0) {
  pcol <- prob_col_name(outer_pred_dt)
  if (pcol == "prob.case") {
    outer_pred_dt[, prob_case := `prob.case`]
  }

  fwrite(
    outer_pred_dt,
    file.path(outercv_dir, paste0("outercv_prediction_level_", disease_id, ".csv"))
  )

  outercv_metrics_ci <- outer_pred_dt[
    ,
    calc_metrics_by_config(.SD),
    by = .(MLA, HC, FSF, CBA)
  ]

  outercv_metrics_pretty <- outercv_metrics_ci |>
    as_tibble() |>
    mutate(
      AUROC_CI = fmt_ci(AUROC, AUROC_l, AUROC_u),
      BACC_CI  = fmt_ci(BACC, BACC_l, BACC_u),
      Sens_CI  = fmt_ci(Sens, Sens_l, Sens_u),
      Spec_CI  = fmt_ci(Spec, Spec_l, Spec_u),
      Brier_CI = fmt_ci(Brier, Brier_l, Brier_u),
      Configuration = paste(MLA, HC, FSF, CBA, sep = " | ")
    ) |>
    arrange(desc(AUROC), desc(BACC))

  fwrite(
    as.data.table(outercv_metrics_ci),
    file.path(outercv_dir, "outercv_metrics_with_ci.csv")
  )

  fwrite(
    as.data.table(outercv_metrics_pretty),
    file.path(outercv_dir, "outercv_metrics_with_ci_pretty.csv")
  )
}

# ==============================================================================
# 13. Permutation test for selected configurations
# ==============================================================================

run_perm_test_config <- function(df, n_perm = 1000, seed = 7832) {

  set.seed(seed)

  config_name <- unique(as.character(df$Configuration))

  roc_obs <- pROC::roc(
    response = df$truth,
    predictor = df$prob_case,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )

  auc_obs <- as.numeric(pROC::auc(roc_obs))

  auc_null <- replicate(n_perm, {
    perm_truth <- sample(df$truth)

    as.numeric(
      pROC::auc(
        pROC::roc(
          response = perm_truth,
          predictor = df$prob_case,
          levels = c(negative_class, positive_class),
          direction = "<",
          quiet = TRUE
        )
      )
    )
  })

  p_val <- (sum(auc_null >= auc_obs) + 1) / (n_perm + 1)

  tibble(
    Configuration = config_name,
    Observed_AUROC = auc_obs,
    Null_Mean_AUROC = mean(auc_null),
    Null_SD_AUROC = sd(auc_null),
    Null_95CI = sprintf(
      "%.3f--%.3f",
      quantile(auc_null, 0.025),
      quantile(auc_null, 0.975)
    ),
    Permutation_p = p_val,
    n_permutations = n_perm
  )
}

perm_source <- if (exists("outer_pred_dt") && nrow(outer_pred_dt) > 0) outer_pred_dt else pred_dt

perm_dt <- perm_source |>
  as_tibble() |>
  filter(Configuration %in% top4_configs_for_permutation) |>
  mutate(Configuration = factor(Configuration, levels = top4_configs_for_permutation))

if (nrow(perm_dt) > 0) {
  perm_table <- perm_dt |>
    group_by(Configuration) |>
    group_split() |>
    map_dfr(~ run_perm_test_config(.x, n_perm = 1000, seed = 7832)) |>
    mutate(
      Configuration = factor(Configuration, levels = top4_configs_for_permutation)
    ) |>
    arrange(Configuration) |>
    mutate(
      Observed_AUROC = sprintf("%.3f", Observed_AUROC),
      Null_Mean_AUROC = sprintf("%.3f", Null_Mean_AUROC),
      Null_SD_AUROC = sprintf("%.3f", Null_SD_AUROC),
      Permutation_p = ifelse(
        Permutation_p < 0.001,
        "<0.001",
        sprintf("%.3f", Permutation_p)
      )
    )

  fwrite(
    as.data.table(perm_table),
    file.path(outercv_dir, paste0(disease_id, "_permutation_test_top4_table.csv"))
  )
}

# ==============================================================================
# 14. Cross-cohort validation using tuned learners
# ==============================================================================

ccv_prediction_rows <- list()
failed_runs <- list()

for (f in nestedcv_files) {

  rr <- readRDS(f)
  fname <- basename(f)

  clean_name <- fname |>
    str_remove("^nestedcv_") |>
    str_remove("\\.rds$")

  learner_key <- str_extract(clean_name, "^(RandomForest|LightGBM|RBF_SVM)")
  task_name <- clean_name |>
    str_remove("^(RandomForest|LightGBM|RBF_SVM)_task_")

  model_label <- learners_param_space[[learner_key]]$label
  task_path <- file.path(balanced_task_dir, paste0("train_task_", task_name, ".rds"))

  if (!file.exists(task_path)) {
    failed_runs[[length(failed_runs) + 1]] <- data.table(
      file = fname,
      MLA = model_label,
      pipeline = task_name,
      cohort = NA_character_,
      error = "Task file not found"
    )
    next
  }

  task_base <- readRDS(task_path)

  cohort_dt <- tryCatch(
    as.data.table(task_base$data(cols = "cohort")),
    error = function(e) NULL
  )

  if (is.null(cohort_dt) || !"cohort" %in% names(cohort_dt)) {
    failed_runs[[length(failed_runs) + 1]] <- data.table(
      file = fname,
      MLA = model_label,
      pipeline = task_name,
      cohort = NA_character_,
      error = "cohort could not be extracted"
    )
    next
  }

  cohort_vec <- as.character(cohort_dt$cohort)

  task_base$col_roles$stratum <- character(0)
  task_base$col_roles$group <- character(0)
  task_base$col_roles$name <- character(0)

  if ("cohort" %in% task_base$feature_names) {
    task_base$set_col_roles("cohort", remove_from = "feature")
  }

  best_params <- extract_best_params(rr)
  base_learner <- build_fixed_learner(model_label, best_params)

  cohorts <- sort(unique(na.omit(cohort_vec)))

  for (cohort_i in cohorts) {

    test_ids <- task_base$row_ids[cohort_vec == cohort_i]
    train_ids <- setdiff(task_base$row_ids, test_ids)

    train_y <- task_base$truth(rows = train_ids)

    if (length(unique(train_y)) < 2) {
      failed_runs[[length(failed_runs) + 1]] <- data.table(
        file = fname,
        MLA = model_label,
        pipeline = task_name,
        cohort = cohort_i,
        error = "Single-class training set"
      )
      next
    }

    current_learner <- base_learner$clone(deep = TRUE)

    pred <- tryCatch({
      set.seed(7832)
      current_learner$train(task_base, row_ids = train_ids)
      current_learner$predict(task_base, row_ids = test_ids)
    }, error = function(e) {
      failed_runs[[length(failed_runs) + 1]] <<- data.table(
        file = fname,
        MLA = model_label,
        pipeline = task_name,
        cohort = cohort_i,
        error = e$message
      )
      NULL
    })

    if (!is.null(pred)) {
      dt_row <- as.data.table(pred)
      meta <- parse_task_id(task_name)

      dt_row[, `:=`(
        cohort = cohort_i,
        MLA = model_label,
        pipeline = task_name,
        file = fname,
        HC = meta$HC,
        FSF = meta$FSF,
        CBA = meta$CBA,
        Configuration = paste(model_label, meta$HC, meta$FSF, meta$CBA, sep = " | ")
      )]

      ccv_prediction_rows[[length(ccv_prediction_rows) + 1]] <- dt_row
    }
  }
}

flat_dt <- rbindlist(ccv_prediction_rows, fill = TRUE)

if (nrow(flat_dt) > 0) {
  pcol <- prob_col_name(flat_dt)
  if (pcol == "prob.case") {
    flat_dt[, prob_case := `prob.case`]
  }

  fwrite(
    flat_dt,
    file.path(ccv_dir, paste0(disease_id, "_ccv_prediction_level.csv"))
  )

  saveRDS(
    flat_dt,
    file.path(ccv_dir, paste0(disease_id, "_flat_dt_clean.rds"))
  )
}

if (length(failed_runs) > 0) {
  failed_dt <- rbindlist(failed_runs, fill = TRUE)

  fwrite(
    failed_dt,
    file.path(ccv_dir, paste0(disease_id, "_ccv_failed_runs.csv"))
  )
}

# ==============================================================================
# 15. CCV ROC, PR, calibration and permutation testing for selected configurations
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(pROC)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(readr)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 15.1 User-defined selected CCV configurations
# ------------------------------------------------------------------------------

ccv_fig_dir <- file.path(ccv_dir, "selected_configuration_figures")
ccv_roc_dir <- file.path(ccv_fig_dir, "roc_plots")
ccv_pr_dir  <- file.path(ccv_fig_dir, "pr_plots")
ccv_cal_dir <- file.path(ccv_fig_dir, "calibration_plots")

dir.create(ccv_roc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ccv_pr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ccv_cal_dir, recursive = TRUE, showWarnings = FALSE)

# Edit this table for TB / UC / RA.
# These are the exact configurations plotted in the CCV ROC, PR and calibration figures.
config_levels <- c(
  "RBF SVM | VST + SVA | Info_Gain | SMOTE",
  "Random Forest | VST + SVA | Info_Gain | SMOTE",
  "LightGBM | VST + SVA | Info_Gain | SMOTE",
  "Random Forest | VST + SVA | Hybrid A | SMOTE"
)

# Edit this for your cohort names and desired order.
cohort_levels <- c("Shchetynsky", "Guo", "Sales")

curve_pal <- c(
  "Shchetynsky" = "#E41A1C",
  "Guo" = "#377EB8",
  "Sales" = "#4DAF4A",
  "All cohorts" = "black"
)

# ------------------------------------------------------------------------------
# 15.2 Standardise CCV prediction table
# ------------------------------------------------------------------------------

if (!exists("flat_dt") || nrow(flat_dt) == 0) {
  stop("flat_dt was not found or is empty. Run the CCV loop before Step 15.", call. = FALSE)
}

setDT(flat_dt)

if (!"prob.case" %in% names(flat_dt) && "prob_case" %in% names(flat_dt)) {
  flat_dt[, `prob.case` := prob_case]
}

required_ccv_cols <- c("truth", "prob.case", "cohort", "MLA", "HC", "FSF", "CBA")
missing_ccv_cols <- setdiff(required_ccv_cols, names(flat_dt))

if (length(missing_ccv_cols) > 0) {
  stop(
    "flat_dt is missing required column(s): ",
    paste(missing_ccv_cols, collapse = ", "),
    call. = FALSE
  )
}

all_preds_final <- flat_dt %>%
  as_tibble() %>%
  mutate(
    MLA = case_when(
      str_detect(as.character(file), "RandomForest") ~ "Random Forest",
      str_detect(as.character(file), "RBF_SVM") ~ "RBF SVM",
      str_detect(as.character(file), "LightGBM") ~ "LightGBM",
      TRUE ~ as.character(MLA)
    ),

    HC = case_when(
      str_detect(as.character(file), "VST_ComBatRef") ~ "VST + ComBat-Ref",
      str_detect(as.character(file), "VST_ComBat") ~ "VST + ComBat",
      str_detect(as.character(file), "VST_SVA") ~ "VST + SVA",
      str_detect(as.character(file), "VST_only") ~ "VST only",
      str_detect(as.character(file), "ComBatSeq_log1p") ~ "ComBat-Seq + log1p",
      TRUE ~ as.character(HC)
    ),

    FSF = case_when(
      str_detect(as.character(file), "HybridA") ~ "Hybrid A",
      str_detect(as.character(file), "HybridB") ~ "Hybrid B",
      str_detect(as.character(file), "information_gain") ~ "Info_Gain",
      str_detect(as.character(file), "boruta") ~ "Boruta",
      str_detect(as.character(file), "mrmr") ~ "mRMR",
      TRUE ~ as.character(FSF)
    ),

    CBA = case_when(
      str_detect(as.character(file), "SMOTE") ~ "SMOTE",
      str_detect(as.character(file), "Downsample") ~ "Downsample",
      str_detect(as.character(file), "Upsample") ~ "Upsample",
      TRUE ~ as.character(CBA)
    ),

    cohort = str_trim(as.character(cohort)),
    truth = factor(as.character(truth), levels = c(negative_class, positive_class)),
    y = ifelse(as.character(truth) == positive_class, 1, 0),
    p = as.numeric(`prob.case`),

    Configuration = paste(MLA, HC, FSF, CBA, sep = " | ")
  )

all_preds_final %>%
  count(Configuration, sort = TRUE) %>%
  print(n = 50)

subset_pred <- all_preds_final %>%
  filter(
    Configuration %in% config_levels,
    !is.na(p)
  ) %>%
  group_by(Configuration, cohort) %>%
  filter(
    n() >= 10,
    n_distinct(truth) == 2
  ) %>%
  ungroup() %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = cohort_levels)
  ) %>%
  filter(!is.na(Configuration), !is.na(cohort)) %>%
  droplevels()

print(table(subset_pred$Configuration, useNA = "ifany"))
print(table(subset_pred$Configuration, subset_pred$cohort))

if (nrow(subset_pred) == 0) {
  stop("No CCV predictions matched config_levels. Check MLA, HC, FSF and CBA labels.", call. = FALSE)
}

# ==============================================================================
# 15A. ROC curves: cohort-specific + pooled all-cohorts line
# ==============================================================================

create_roc_df <- function(df,
                          cohort_name,
                          fpr_grid = seq(0, 1, by = 0.01),
                          boot_n = 300) {

  roc_obj <- pROC::roc(
    response = df$truth,
    predictor = df$p,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )

  raw_roc <- pROC::coords(
    roc_obj,
    x = "all",
    ret = c("specificity", "sensitivity"),
    transpose = FALSE
  ) %>%
    as_tibble() %>%
    transmute(
      fpr = 1 - as.numeric(specificity),
      sensitivity = as.numeric(sensitivity)
    ) %>%
    bind_rows(
      tibble(fpr = 0, sensitivity = 0),
      tibble(fpr = 1, sensitivity = 1)
    ) %>%
    arrange(fpr, sensitivity) %>%
    distinct(fpr, .keep_all = TRUE) %>%
    mutate(sensitivity = cummax(sensitivity))

  interp_sens <- approx(
    x = raw_roc$fpr,
    y = raw_roc$sensitivity,
    xout = fpr_grid,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y

  roc_df <- tibble(
    fpr = fpr_grid,
    sensitivity = pmin(pmax(interp_sens, 0), 1),
    cohort = cohort_name,
    Configuration = as.character(unique(df$Configuration))
  )

  ci_obj <- try(
    pROC::ci.se(
      roc_obj,
      specificities = 1 - fpr_grid,
      boot.n = boot_n
    ),
    silent = TRUE
  )

  if (!inherits(ci_obj, "try-error")) {
    ci_tmp <- as.data.frame(ci_obj)
    colnames(ci_tmp) <- c("sens_lower", "sensitivity_mid", "sens_upper")

    ci_df <- tibble(
      fpr = fpr_grid,
      sens_lower = pmax(as.numeric(ci_tmp$sens_lower), 0),
      sens_upper = pmin(as.numeric(ci_tmp$sens_upper), 1),
      cohort = cohort_name,
      Configuration = as.character(unique(df$Configuration))
    )
  } else {
    ci_df <- roc_df %>%
      transmute(
        fpr,
        sens_lower = sensitivity,
        sens_upper = sensitivity,
        cohort,
        Configuration
      )
  }

  auc_ci <- try(pROC::ci.auc(roc_obj, method = "delong"), silent = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))

  if (!inherits(auc_ci, "try-error")) {
    auc_low <- as.numeric(auc_ci[1])
    auc_high <- as.numeric(auc_ci[3])
  } else {
    auc_low <- NA_real_
    auc_high <- NA_real_
  }

  auc_df <- tibble(
    Configuration = as.character(unique(df$Configuration)),
    cohort = cohort_name,
    auc = auc_val,
    low = auc_low,
    high = auc_high,
    cohort_label = ifelse(
      is.na(auc_low) | is.na(auc_high),
      sprintf("%s (AUROC %.3f)", cohort_name, auc_val),
      sprintf(
        "%s (AUROC %.3f [%.3f–%.3f])",
        cohort_name, auc_val, auc_low, auc_high
      )
    )
  )

  list(
    roc_df = roc_df,
    ci_df = ci_df,
    auc_df = auc_df
  )
}

roc_results <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_split() %>%
  map(~ create_roc_df(
    df = .x,
    cohort_name = as.character(unique(.x$cohort))
  ))

roc_df_all <- map_dfr(roc_results, "roc_df")
ci_df_all <- map_dfr(roc_results, "ci_df")
auc_df_all <- map_dfr(roc_results, "auc_df")

all_results <- subset_pred %>%
  group_by(Configuration) %>%
  group_split() %>%
  map(~ create_roc_df(
    df = .x,
    cohort_name = "All cohorts"
  ))

all_roc_df <- map_dfr(all_results, "roc_df")
all_ci_df <- map_dfr(all_results, "ci_df")
all_auc_df <- map_dfr(all_results, "auc_df")

roc_table <- bind_rows(
  auc_df_all,
  all_auc_df
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  ) %>%
  arrange(Configuration, cohort)

write_csv(
  roc_table,
  file.path(ccv_roc_dir, "ROC_table_with_all_cohorts_DeLong95CI.csv")
)

roc_table_wide <- roc_table %>%
  select(Configuration, cohort, cohort_label) %>%
  pivot_wider(
    names_from = cohort,
    values_from = cohort_label
  )

write_csv(
  roc_table_wide,
  file.path(ccv_roc_dir, "ROC_table_wide_with_all_cohorts_DeLong95CI.csv")
)

auc_panel_labels <- roc_table %>%
  mutate(
    cohort = as.character(cohort),
    label = sprintf(
      "%s: %.3f [%.3f–%.3f]",
      cohort, auc, low, high
    ),
    y = case_when(
      cohort == cohort_levels[1] ~ 0.18,
      cohort == cohort_levels[2] ~ 0.13,
      cohort == cohort_levels[3] ~ 0.08,
      cohort == "All cohorts" ~ 0.03,
      TRUE ~ 0.03
    ),
    x_line_start = 0.40,
    x_line_end = 0.47,
    x_text = 0.49
  ) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  )

roc_plot_df <- bind_rows(
  roc_df_all,
  all_roc_df
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  )

ci_plot_df <- bind_rows(
  ci_df_all,
  all_ci_df
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  )

p_roc_all <- ggplot() +
  geom_ribbon(
    data = ci_plot_df %>% filter(cohort != "All cohorts"),
    aes(
      x = fpr,
      ymin = sens_lower,
      ymax = sens_upper,
      fill = cohort
    ),
    alpha = 0.14,
    colour = NA
  ) +
  geom_ribbon(
    data = ci_plot_df %>% filter(cohort == "All cohorts"),
    aes(
      x = fpr,
      ymin = sens_lower,
      ymax = sens_upper
    ),
    fill = "grey40",
    alpha = 0.16,
    colour = NA
  ) +
  geom_line(
    data = roc_plot_df %>% filter(cohort != "All cohorts"),
    aes(
      x = fpr,
      y = sensitivity,
      colour = cohort,
      group = interaction(Configuration, cohort)
    ),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_line(
    data = roc_plot_df %>% filter(cohort == "All cohorts"),
    aes(
      x = fpr,
      y = sensitivity,
      colour = cohort,
      group = interaction(Configuration, cohort)
    ),
    linewidth = 1.35,
    lineend = "round"
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.7
  ) +
  geom_rect(
    data = auc_panel_labels %>%
      distinct(Configuration),
    aes(
      xmin = 0.36,
      xmax = 1.02,
      ymin = -0.01,
      ymax = 0.24
    ),
    inherit.aes = FALSE,
    fill = alpha("white", 0.92),
    colour = "black",
    linewidth = 0.45
  ) +
  geom_segment(
    data = auc_panel_labels,
    aes(
      x = x_line_start,
      xend = x_line_end,
      y = y,
      yend = y,
      colour = cohort
    ),
    inherit.aes = FALSE,
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_text(
    data = auc_panel_labels,
    aes(
      x = x_text,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.45,
    colour = "black",
    lineheight = 0.95
  ) +
  facet_wrap(~ Configuration, nrow = 2, ncol = 2, drop = FALSE) +
  scale_colour_manual(
    values = curve_pal,
    breaks = c(cohort_levels, "All cohorts"),
    name = "ROC curve"
  ) +
  scale_fill_manual(
    values = curve_pal,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  labs(
    title = paste0("Cross-cohort ROC curves for selected ", disease_id, " configurations"),
    subtitle = "Coloured lines show held-out cohorts; black line shows pooled held-out predictions",
    x = "False Positive Rate (1 − Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 8, margin = margin(t = 4)),
    axis.text.y = element_text(colour = "black", size = 8, margin = margin(r = 4)),
    axis.title = element_text(colour = "black", face = "bold", size = 10.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5),
    plot.subtitle = element_text(colour = "black", size = 10.5, hjust = 0.5),
    strip.text = element_text(face = "bold", colour = "black", size = 8.2),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 16, 8, 16)
  )

print(p_roc_all)

ggsave(
  filename = file.path(ccv_roc_dir, paste0(disease_id, "_CCV_ROC_selected_configs.png")),
  plot = p_roc_all,
  width = 8.5,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(ccv_roc_dir, paste0(disease_id, "_CCV_ROC_selected_configs.pdf")),
  plot = p_roc_all,
  width = 8.5,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

print(roc_table)

# ==============================================================================
# 15B. PR curves: cohort-specific + pooled all-cohorts line
# ==============================================================================

auprc_once <- function(y, p) {
  o <- order(-p)
  y <- y[o]

  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)

  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / sum(y == 1)

  recall <- c(0, recall)
  precision <- c(1, precision)

  sum(diff(recall) * precision[-length(precision)], na.rm = TRUE)
}

boot_auprc <- function(df, B = 1000, seed = 7832) {
  set.seed(seed)

  n <- nrow(df)
  est <- auprc_once(df$y, df$p)

  boots <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    auprc_once(df$y[idx], df$p[idx])
  })

  ci <- quantile(boots, c(0.025, 0.975), names = FALSE, na.rm = TRUE)

  c(est = est, low = ci[1], high = ci[2])
}

pr_points_one <- function(df) {
  o <- order(-df$p)
  y <- df$y[o]

  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)

  tibble(
    recall = c(0, tp / sum(y == 1)),
    precision = c(1, tp / pmax(tp + fp, 1))
  ) %>%
    arrange(recall, desc(precision)) %>%
    distinct(recall, .keep_all = TRUE)
}

boot_pr_band_one <- function(df,
                             B = 300,
                             seed = 7832,
                             recall_grid = seq(0, 1, by = 0.01)) {
  set.seed(seed)

  n <- nrow(df)

  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    prb <- pr_points_one(df[idx, , drop = FALSE])

    approx(
      x = prb$recall,
      y = prb$precision,
      xout = recall_grid,
      method = "linear",
      ties = "ordered",
      rule = 2
    )$y
  })

  tibble(
    recall = recall_grid,
    precision_low = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    precision_mid = apply(boot_mat, 1, median, na.rm = TRUE),
    precision_high = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

pr_points <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_modify(~ pr_points_one(.x)) %>%
  ungroup()

pr_band <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_modify(~ boot_pr_band_one(.x, B = 300, seed = 7832)) %>%
  ungroup()

auprc_tbl <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  summarise(
    prevalence = mean(y),
    tmp = list(boot_auprc(cur_data(), B = 1000, seed = 7832)),
    .groups = "drop"
  ) %>%
  mutate(
    AUPRC = map_dbl(tmp, 1),
    low = map_dbl(tmp, 2),
    high = map_dbl(tmp, 3),
    cohort_label = sprintf(
      "%s (AUPRC %.3f [%.3f–%.3f])",
      cohort, AUPRC, low, high
    )
  ) %>%
  select(-tmp)

all_pr_points <- subset_pred %>%
  group_by(Configuration) %>%
  group_modify(~ pr_points_one(.x)) %>%
  ungroup() %>%
  mutate(cohort = "All cohorts")

all_pr_band <- subset_pred %>%
  group_by(Configuration) %>%
  group_modify(~ boot_pr_band_one(.x, B = 300, seed = 7832)) %>%
  ungroup() %>%
  mutate(cohort = "All cohorts")

all_auprc_tbl <- subset_pred %>%
  group_by(Configuration) %>%
  summarise(
    prevalence = mean(y),
    tmp = list(boot_auprc(cur_data(), B = 1000, seed = 7832)),
    .groups = "drop"
  ) %>%
  mutate(
    cohort = "All cohorts",
    AUPRC = map_dbl(tmp, 1),
    low = map_dbl(tmp, 2),
    high = map_dbl(tmp, 3),
    cohort_label = sprintf(
      "%s (AUPRC %.3f [%.3f–%.3f])",
      cohort, AUPRC, low, high
    )
  ) %>%
  select(Configuration, cohort, prevalence, AUPRC, low, high, cohort_label)

pr_table <- bind_rows(
  auprc_tbl,
  all_auprc_tbl
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  ) %>%
  arrange(Configuration, cohort)

write_csv(
  pr_table,
  file.path(ccv_pr_dir, "PR_table_with_all_cohorts_bootstrap95CI.csv")
)

pr_table_wide <- pr_table %>%
  select(Configuration, cohort, cohort_label) %>%
  pivot_wider(
    names_from = cohort,
    values_from = cohort_label
  )

write_csv(
  pr_table_wide,
  file.path(ccv_pr_dir, "PR_table_wide_with_all_cohorts_bootstrap95CI.csv")
)

pr_panel_labels <- pr_table %>%
  mutate(
    cohort = as.character(cohort),
    label = sprintf(
      "%s: %.3f [%.3f–%.3f]",
      cohort, AUPRC, low, high
    ),
    y = case_when(
      cohort == cohort_levels[1] ~ 0.14,
      cohort == cohort_levels[2] ~ 0.10,
      cohort == cohort_levels[3] ~ 0.06,
      cohort == "All cohorts" ~ 0.02,
      TRUE ~ 0.02
    ),
    x_line_start = 0.48,
    x_line_end = 0.55,
    x_text = 0.57
  ) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  ) %>%
  filter(!is.na(Configuration))

pr_points_all <- bind_rows(
  pr_points,
  all_pr_points
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  )

pr_band_all <- bind_rows(
  pr_band,
  all_pr_band
) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels),
    cohort = factor(cohort, levels = c(cohort_levels, "All cohorts"))
  )

p_pr_all <- ggplot() +
  geom_ribbon(
    data = pr_band_all %>% filter(cohort != "All cohorts"),
    aes(
      x = recall,
      ymin = precision_low,
      ymax = precision_high,
      fill = cohort
    ),
    alpha = 0.14,
    colour = NA
  ) +
  geom_ribbon(
    data = pr_band_all %>% filter(cohort == "All cohorts"),
    aes(
      x = recall,
      ymin = precision_low,
      ymax = precision_high
    ),
    fill = "grey40",
    alpha = 0.16,
    colour = NA
  ) +
  geom_line(
    data = pr_points_all %>% filter(cohort != "All cohorts"),
    aes(
      x = recall,
      y = precision,
      colour = cohort,
      group = interaction(Configuration, cohort)
    ),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_line(
    data = pr_points_all %>% filter(cohort == "All cohorts"),
    aes(
      x = recall,
      y = precision,
      colour = cohort,
      group = interaction(Configuration, cohort)
    ),
    linewidth = 1.35,
    lineend = "round"
  ) +
  geom_hline(
    data = auprc_tbl,
    aes(
      yintercept = prevalence,
      colour = cohort
    ),
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_hline(
    data = all_auprc_tbl,
    aes(yintercept = prevalence),
    colour = "black",
    linetype = "dashed",
    linewidth = 0.9,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_segment(
    data = pr_panel_labels,
    aes(
      x = x_line_start,
      xend = x_line_end,
      y = y,
      yend = y,
      colour = cohort
    ),
    inherit.aes = FALSE,
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_text(
    data = pr_panel_labels,
    aes(
      x = x_text,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.0,
    colour = "black",
    lineheight = 0.9
  ) +
  facet_wrap(~ Configuration, nrow = 2, ncol = 2, drop = FALSE) +
  scale_colour_manual(
    values = curve_pal,
    breaks = c(cohort_levels, "All cohorts"),
    name = "PR curve"
  ) +
  scale_fill_manual(
    values = curve_pal,
    guide = "none"
  ) +
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
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  labs(
    title = paste0("Cross-cohort precision-recall curves for selected ", disease_id, " configurations"),
    subtitle = "Coloured lines show held-out cohorts; black line shows pooled held-out predictions",
    x = "Recall (Sensitivity)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 8, margin = margin(t = 5)),
    axis.text.y = element_text(colour = "black", size = 8, margin = margin(r = 5)),
    axis.title = element_text(colour = "black", face = "bold", size = 10.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    plot.title = element_text(face = "bold", colour = "black", size = 15, hjust = 0.5),
    plot.subtitle = element_text(colour = "black", size = 10.5, hjust = 0.5),
    strip.text = element_text(face = "bold", colour = "black", size = 9.5),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 22, 8, 22)
  )

print(p_pr_all)

ggsave(
  filename = file.path(ccv_pr_dir, paste0(disease_id, "_CCV_PR_selected_configs.png")),
  plot = p_pr_all,
  width = 8.5,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(ccv_pr_dir, paste0(disease_id, "_CCV_PR_selected_configs.pdf")),
  plot = p_pr_all,
  width = 8.5,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

print(pr_table)

# ==============================================================================
# 15C. CCV calibration curves
# ==============================================================================

calc_calibration <- function(df, n_bins = 10) {

  d <- df %>%
    mutate(bin = ntile(p, n_bins))

  bins <- d %>%
    group_by(bin) %>%
    summarise(
      conf = mean(p),
      acc = mean(y),
      n = n(),
      acc_lo = binom.test(sum(y), n)$conf.int[1],
      acc_hi = binom.test(sum(y), n)$conf.int[2],
      .groups = "drop"
    )

  ece <- sum((bins$n / sum(bins$n)) * abs(bins$acc - bins$conf))
  mce <- max(abs(bins$acc - bins$conf))
  brier <- mean((d$p - d$y)^2)

  list(
    bins = bins,
    ECE = ece,
    MCE = mce,
    Brier = brier
  )
}

calibration_bins <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_modify(~ {
    tmp <- calc_calibration(.x)

    tmp$bins %>%
      mutate(
        ECE = tmp$ECE,
        MCE = tmp$MCE,
        Brier = tmp$Brier
      )
  }) %>%
  ungroup()

loess_lines <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_modify(~ {
    fit <- loess(y ~ p, data = .x, span = 0.75)

    grid <- tibble(
      p = seq(0, 1, by = 0.01)
    )

    grid$m <- predict(fit, newdata = grid)

    grid %>%
      mutate(
        m = pmin(pmax(m, 0), 1)
      )
  }) %>%
  ungroup()

boot_calibration_band <- function(df,
                                  B = 300,
                                  grid = seq(0, 1, by = 0.01),
                                  span = 0.75,
                                  seed = 7832) {

  set.seed(seed)
  n <- nrow(df)

  mats <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    d_b <- df[idx, ]

    fit <- try(
      loess(y ~ p, data = d_b, span = span),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(rep(NA_real_, length(grid)))
    }

    predict(
      fit,
      newdata = data.frame(p = grid)
    )
  })

  tibble(
    p = grid,
    lo = apply(mats, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi = apply(mats, 1, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      lo = pmax(pmin(lo, 1), 0),
      hi = pmax(pmin(hi, 1), 0)
    )
}

calibration_bands <- subset_pred %>%
  group_by(Configuration, cohort) %>%
  group_modify(~ boot_calibration_band(.x)) %>%
  ungroup()

cal_labels <- calibration_bins %>%
  group_by(Configuration, cohort) %>%
  summarise(
    Brier = first(Brier),
    ECE = first(ECE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "%s\nBrier %.3f; ECE %.3f",
      cohort,
      Brier,
      ECE
    ),
    y = case_when(
      as.character(cohort) == cohort_levels[1] ~ 0.18,
      as.character(cohort) == cohort_levels[2] ~ 0.11,
      as.character(cohort) == cohort_levels[3] ~ 0.04,
      TRUE ~ 0.04
    ),
    x_line_start = 0.62,
    x_line_end = 0.70,
    x_text = 0.72
  )

p_cal_ccv <- ggplot() +
  geom_ribbon(
    data = calibration_bands,
    aes(
      x = p,
      ymin = lo,
      ymax = hi,
      fill = cohort
    ),
    alpha = 0.14,
    colour = NA
  ) +
  geom_line(
    data = loess_lines,
    aes(
      x = p,
      y = m,
      colour = cohort,
      group = interaction(Configuration, cohort)
    ),
    linewidth = 1.05
  ) +
  geom_errorbar(
    data = calibration_bins,
    aes(
      x = conf,
      ymin = acc_lo,
      ymax = acc_hi,
      colour = cohort
    ),
    width = 0.02,
    linewidth = 0.4
  ) +
  geom_point(
    data = calibration_bins,
    aes(
      x = conf,
      y = acc,
      colour = cohort
    ),
    size = 1.8
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_segment(
    data = cal_labels,
    aes(
      x = x_line_start,
      xend = x_line_end,
      y = y,
      yend = y,
      colour = cohort
    ),
    linewidth = 1.0
  ) +
  geom_text(
    data = cal_labels,
    aes(
      x = x_text,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 2.0,
    lineheight = 0.9
  ) +
  facet_wrap(
    ~ Configuration,
    nrow = 2,
    ncol = 2,
    drop = FALSE
  ) +
  scale_colour_manual(values = curve_pal) +
  scale_fill_manual(values = curve_pal) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = percent_format(accuracy = 1)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  labs(
    title = paste0("Cross-cohort calibration curves for selected ", disease_id, " configurations"),
    subtitle = "Held-out cohort calibration under independent cohort validation",
    x = "Predicted probability of case",
    y = "Observed fraction of cases"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(colour = "black", size = 8),
    axis.title = element_text(face = "bold", size = 10.5),
    strip.text = element_text(face = "bold", size = 9.5),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5)
  )

print(p_cal_ccv)

ggsave(
  filename = file.path(ccv_cal_dir, paste0(disease_id, "_CCV_calibration_selected_configs.png")),
  plot = p_cal_ccv,
  width = 8.5,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(ccv_cal_dir, paste0(disease_id, "_CCV_calibration_selected_configs.pdf")),
  plot = p_cal_ccv,
  width = 8.5,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

# ==============================================================================
# 15D. Permutation test for selected CCV configurations
# ==============================================================================

perm_dt <- subset_pred %>%
  filter(Configuration %in% config_levels) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels)
  )

run_perm_test_config <- function(df, n_perm = 1000, seed = 7832) {

  set.seed(seed)

  config_name <- unique(as.character(df$Configuration))

  roc_obs <- pROC::roc(
    response = df$truth,
    predictor = df$p,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )

  auc_obs <- as.numeric(pROC::auc(roc_obs))

  auc_null <- replicate(n_perm, {
    perm_truth <- sample(df$truth)

    as.numeric(
      pROC::auc(
        pROC::roc(
          response = perm_truth,
          predictor = df$p,
          levels = c(negative_class, positive_class),
          direction = "<",
          quiet = TRUE
        )
      )
    )
  })

  p_val <- (sum(auc_null >= auc_obs) + 1) / (n_perm + 1)

  tibble(
    Configuration = config_name,
    Observed_AUROC = auc_obs,
    Null_Mean_AUROC = mean(auc_null),
    Null_SD_AUROC = sd(auc_null),
    Null_95CI = sprintf(
      "%.3f--%.3f",
      quantile(auc_null, 0.025),
      quantile(auc_null, 0.975)
    ),
    Permutation_p = p_val,
    n_permutations = n_perm
  )
}

perm_table <- perm_dt %>%
  group_by(Configuration) %>%
  group_split() %>%
  map_dfr(~ run_perm_test_config(.x, n_perm = 1000, seed = 7832)) %>%
  mutate(
    Configuration = factor(Configuration, levels = config_levels)
  ) %>%
  arrange(Configuration) %>%
  mutate(
    Observed_AUROC = sprintf("%.3f", Observed_AUROC),
    Null_Mean_AUROC = sprintf("%.3f", Null_Mean_AUROC),
    Null_SD_AUROC = sprintf("%.3f", Null_SD_AUROC),
    Permutation_p = ifelse(
      as.numeric(Permutation_p) < 0.001,
      "<0.001",
      sprintf("%.3f", as.numeric(Permutation_p))
    )
  )

print(perm_table)

fwrite(
  as.data.table(perm_table),
  file.path(ccv_fig_dir, paste0(disease_id, "_CCV_permutation_test_config_levels_table.csv"))
)

# ==============================================================================
# 16. Save reproducibility log
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(summary_dir, paste0(disease_id, "_sessionInfo.txt"))
)

message("\n[done] Script 03 completed successfully.")
message("[done] Baseline outputs: ", baseline_dir)
message("[done] Nested CV outputs: ", ncv_dir)
message("[done] Outer CV outputs: ", outercv_dir)
message("[done] CCV outputs: ", ccv_dir)
message("[done] Figures: ", fig_dir)
