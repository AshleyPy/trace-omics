# ==============================================================================
# TRACE-Omics
# Script 02: Feature-selection framework comparisons and feature-space conditioning
#
# Purpose:
#   This script applies single feature-selection filters to pre-saved mlr3 training
#   tasks, compares selected feature sets using Jaccard Index, and quantifies
#   feature-space conditioning using effective rank and variance inflation factors.
#
# Required input structure:
#
#   data/processed/TB/tasks/
#   ├── train_task_TB_VST_only.rds
#   ├── train_task_TB_VST_SVA.rds
#   ├── train_task_TB_VST_ComBat.rds
#   ├── train_task_TB_VST_ComBatRef.rds
#   └── train_task_TB_ComBatSeq_log1p.rds
#
# Optional hybrid input structure:
#
#   data/processed/TB/hybrid_features/
#   ├── hybrid_features_HybridA_TB_VST_only.rds
#   ├── hybrid_features_HybridB_TB_VST_only.rds
#   └── ...
#
# Outputs:
#   - Top-N feature sets per filter and harmonisation configuration
#   - Jaccard Index tables and heatmaps
#   - eRank and VIF stability trade-off tables
#   - sessionInfo log
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
  library(mlr3)
  library(mlr3filters)
  library(praznik)
  library(FSelectorRcpp)
  library(mlr3learners)
  library(car)
})

set.seed(7832)

# ==============================================================================
# 1. User settings
# ==============================================================================

disease_id <- "TB"

task_dir <- file.path("data", "processed", disease_id, "tasks")
feature_dir <- file.path("results", "02_fsf_stability_conditioning", disease_id, "feature_subsets")
fig_dir <- file.path("results", "02_fsf_stability_conditioning", disease_id, "figures")
table_dir <- file.path("results", "02_fsf_stability_conditioning", disease_id, "tables")

hybrid_dir <- file.path("data", "processed", disease_id, "hybrid_features")

dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 50
n_boot <- 200

filter_ids <- c(
  "anova",
  "auc",
  "importance",
  "boruta",
  "cmim",
  "mrmr",
  "information_gain",
  "variance"
)

hc_order <- c(
  "VST_only",
  "VST_SVA",
  "VST_ComBat",
  "VST_ComBatRef",
  "ComBatSeq_log1p"
)

name_map <- c(
  anova = "ANOVA",
  auc = "AUC",
  information_gain = "Info_Gain",
  boruta = "Boruta",
  cmim = "CMIM",
  mrmr = "mRMR",
  variance = "RF_Var",
  importance = "RF_Imp"
)

method_cols <- c(
  "ANOVA" = "#E41A1C",
  "AUC" = "#377EB8",
  "Boruta" = "#4DAF4A",
  "CMIM" = "#984EA3",
  "Info_Gain" = "#FF7F00",
  "mRMR" = "#FFFF33",
  "RF_Imp" = "#A65628",
  "RF_Var" = "#F781BF",
  "HybridA" = "steelblue",
  "HybridB" = "firebrick"
)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

safe_read_vec <- function(path) {
  if (!file.exists(path)) return(character(0))

  x <- readRDS(path)

  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }

  if (is.data.frame(x)) {
    x <- x[[1]]
  }

  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}

clean_task_id <- function(path) {
  gsub("^train_task_|\\.rds$", "", basename(path))
}

clean_hc_label <- function(task_id, disease_id) {
  x <- gsub(paste0("^", disease_id, "_"), "", task_id)
  x <- gsub(paste0("^ML_", disease_id, "_"), "", x)
  x
}

clean_filter_key <- function(x) {
  gsub("\\.rds$", "", gsub("^features_", "", x))
}

ji_matrix <- function(sets) {
  nm <- names(sets)
  n <- length(sets)

  M <- matrix(0, n, n, dimnames = list(nm, nm))

  for (i in seq_len(n)) {
    ai <- unique(sets[[i]])

    for (j in seq_len(n)) {
      aj <- unique(sets[[j]])
      U <- union(ai, aj)

      if (length(U) == 0L) {
        M[i, j] <- NA_real_
      } else {
        M[i, j] <- length(intersect(ai, aj)) / length(U)
      }
    }
  }

  diag(M) <- 1
  M
}

order_by_cluster <- function(M) {
  M[is.na(M)] <- 0
  M <- pmin(pmax(M, 0), 1)

  hc <- hclust(as.dist(1 - M), method = "average")
  ord <- hc$order

  M[ord, ord, drop = FALSE]
}

load_feature_sets_for_task <- function(task_feature_dir) {
  files <- list.files(task_feature_dir, "^features_(.+)\\.rds$", full.names = TRUE)

  keep_pat <- paste0("^features_(", paste(filter_ids, collapse = "|"), ")\\.rds$")
  files <- files[grepl(keep_pat, basename(files))]

  if (!length(files)) return(NULL)

  sets <- lapply(files, safe_read_vec)
  names(sets) <- clean_filter_key(basename(files))

  sets[lengths(sets) > 0]
}

ji_tidy_for_task <- function(task_id) {
  task_feature_dir <- file.path(feature_dir, task_id)
  sets <- load_feature_sets_for_task(task_feature_dir)

  if (is.null(sets) || length(sets) < 2L) return(NULL)

  M <- ji_matrix(sets)
  M <- order_by_cluster(M)

  lab <- matrix("", nrow(M), ncol(M), dimnames = dimnames(M))
  lab[lower.tri(M, diag = TRUE)] <- sprintf("%.2f", M[lower.tri(M, diag = TRUE)])

  df_num <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  colnames(df_num) <- c("Filter1", "Filter2", "JI")

  df_lab <- as.data.frame(as.table(lab), stringsAsFactors = FALSE)
  colnames(df_lab) <- c("Filter1", "Filter2", "label")

  left_join(df_num, df_lab, by = c("Filter1", "Filter2")) |>
    mutate(
      task_id = task_id,
      HC = clean_hc_label(task_id, disease_id),
      Filter1 = factor(Filter1, levels = rownames(M)),
      Filter2 = factor(Filter2, levels = colnames(M))
    ) |>
    as_tibble()
}

prep_X <- function(data, feats) {
  feats <- intersect(as.character(feats), colnames(data))

  if (length(feats) < 2L) return(NULL)

  X <- data[, feats, drop = FALSE]

  X <- dplyr::mutate(X, dplyr::across(where(is.logical), as.numeric))
  X <- X[, vapply(X, function(z) is.numeric(z) || is.integer(z), logical(1)), drop = FALSE]

  if (ncol(X) < 2L) return(NULL)

  keep <- vapply(X, function(z) {
    z <- z[is.finite(z)]
    length(z) > 1 && stats::sd(z) > 0
  }, logical(1))

  X <- X[, keep, drop = FALSE]

  if (ncol(X) < 2L) return(NULL)

  mm <- stats::model.matrix(~ 0 + ., data = X)
  qrX <- qr(mm)
  mm <- mm[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]

  mm <- scale(mm)
  attr(mm, "scaled:center") <- NULL
  attr(mm, "scaled:scale") <- NULL

  mm
}

effective_rank <- function(X) {
  s <- svd(X, nu = 0, nv = 0)$d

  if (length(s) == 0) return(NA_real_)

  p <- s^2 / sum(s^2)
  H <- -sum(p * log(p))

  exp(H)
}

mean_vif <- function(X) {
  tryCatch({
    dfX <- as.data.frame(X)

    if (ncol(dfX) < 2L) return(NA_real_)

    mod <- lm(dfX[[1]] ~ ., data = dfX)
    vif_vals <- car::vif(mod)

    mean(vif_vals, na.rm = TRUE)
  }, error = function(e) NA_real_)
}

max_vif <- function(X) {
  tryCatch({
    dfX <- as.data.frame(X)

    if (ncol(dfX) < 2L) return(NA_real_)

    mod <- lm(dfX[[1]] ~ ., data = dfX)
    vif_vals <- car::vif(mod)

    max(vif_vals, na.rm = TRUE)
  }, error = function(e) NA_real_)
}

mean_abs_corr <- function(X) {
  tryCatch({
    cors <- cor(X, use = "pairwise.complete.obs")
    mean(abs(cors[lower.tri(cors)]), na.rm = TRUE)
  }, error = function(e) NA_real_)
}

feature_space_metrics <- function(sets, data, n_boot = 200) {
  results <- list()

  for (nm in names(sets)) {
    feats <- sets[[nm]]
    X <- prep_X(data, feats)

    if (is.null(X)) next

    boot_jacc <- replicate(n_boot, {
      idx <- sample(seq_len(nrow(data)), replace = TRUE)
      subX <- prep_X(data[idx, , drop = FALSE], feats)

      if (is.null(subX)) return(NA_real_)

      intersect_feats <- intersect(colnames(X), colnames(subX))
      union_feats <- union(colnames(X), colnames(subX))

      length(intersect_feats) / length(union_feats)
    })

    boot_erank <- replicate(n_boot, {
      idx <- sample(seq_len(nrow(data)), replace = TRUE)
      subX <- prep_X(data[idx, , drop = FALSE], feats)

      if (is.null(subX)) return(NA_real_)

      effective_rank(subX)
    })

    results[[nm]] <- tibble(
      Filter = nm,
      n_features = length(unique(feats)),
      n_used = ncol(X),
      mean_vif = mean_vif(X),
      max_vif = max_vif(X),
      MeanAbsCorr = mean_abs_corr(X),
      Jaccard = mean(boot_jacc, na.rm = TRUE),
      eRank = mean(boot_erank, na.rm = TRUE),
      eRank_lo = quantile(boot_erank, probs = 0.025, na.rm = TRUE, names = FALSE),
      eRank_hi = quantile(boot_erank, probs = 0.975, na.rm = TRUE, names = FALSE)
    )
  }

  bind_rows(results)
}

load_task_dataframe <- function(task_path) {
  task <- readRDS(task_path)

  df <- as.data.frame(task$data())
  df$outcome <- df[[task$target_names[[1]]]]

  df
}

# ==============================================================================
# 3. Apply single feature-selection filters
# ==============================================================================

train_paths <- list.files(
  task_dir,
  pattern = "^train_task_.*\\.rds$",
  full.names = TRUE
)

if (length(train_paths) == 0) {
  stop("No training tasks found in: ", task_dir, call. = FALSE)
}

message("[tasks] Found ", length(train_paths), " training tasks.")

for (task_path in train_paths) {
  task <- readRDS(task_path)
  task_id <- clean_task_id(task_path)

  message("\n=== Feature selection for task: ", task_id, " ===")

  task_feature_dir <- file.path(feature_dir, task_id)
  dir.create(task_feature_dir, recursive = TRUE, showWarnings = FALSE)

  for (f_id in filter_ids) {
    message("Applying filter: ", f_id)

    top_features <- tryCatch({
      flt_obj <- flt(f_id)
      flt_obj$calculate(task)

      scores_df <- data.frame(
        feature = names(flt_obj$scores),
        score = as.numeric(flt_obj$scores),
        stringsAsFactors = FALSE
      )

      scores_df <- scores_df[order(scores_df$score, decreasing = TRUE), ]

      head(scores_df$feature, top_n)
    }, error = function(e) {
      warning("Filter failed for task ", task_id, " and filter ", f_id, ": ", e$message)
      character(0)
    })

    saveRDS(
      top_features,
      file = file.path(task_feature_dir, paste0("features_", f_id, ".rds"))
    )
  }
}

# ==============================================================================
# 4. Jaccard Index for single filters
# ==============================================================================

message("\n[JI] Computing Jaccard Index matrices.")

task_ids <- clean_task_id(train_paths)

all_ji <- purrr::map_dfr(task_ids, ji_tidy_for_task)

if (nrow(all_ji) > 0) {
  all_ji <- all_ji |>
    mutate(
      HC = factor(HC, levels = hc_order),
      Filter1 = dplyr::recode(as.character(Filter1), !!!name_map),
      Filter2 = dplyr::recode(as.character(Filter2), !!!name_map),
      JI = as.numeric(JI),
      label = ifelse(is.finite(JI), sprintf("%.2f", JI), "")
    ) |>
    filter(Filter1 %in% unname(name_map), Filter2 %in% unname(name_map))

  ord1 <- all_ji |>
    group_by(Filter1) |>
    summarise(m = mean(JI, na.rm = TRUE), .groups = "drop") |>
    rename(method = Filter1)

  ord2 <- all_ji |>
    group_by(Filter2) |>
    summarise(m = mean(JI, na.rm = TRUE), .groups = "drop") |>
    rename(method = Filter2)

  method_order <- bind_rows(ord1, ord2) |>
    group_by(method) |>
    summarise(m = mean(m, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(m)) |>
    pull(method)

  all_ji <- all_ji |>
    mutate(
      Filter1 = factor(Filter1, levels = method_order),
      Filter2 = factor(Filter2, levels = method_order)
    )

  fwrite(
    as.data.table(all_ji),
    file.path(table_dir, paste0(disease_id, "_single_filter_jaccard_index.csv"))
  )

  label_df <- filter(all_ji, nzchar(label))

  p_ji <- ggplot(all_ji, aes(x = Filter2, y = Filter1, fill = JI)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(data = label_df, aes(label = label), size = 2.5, colour = "black") +
    facet_wrap(~ HC, scales = "fixed", ncol = 2) +
    scale_fill_gradientn(
      colours = viridisLite::viridis(256),
      limits = c(0, 1),
      oob = scales::squish,
      name = "Jaccard index"
    ) +
    coord_equal() +
    scale_x_discrete(position = "top") +
    labs(
      title = "Concordance of Feature-Selection Methods",
      subtitle = paste0("Jaccard index of top-", top_n, " selected genes per method"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      legend.position = "right",
      plot.title.position = "plot"
    )

  ggsave(
    file.path(fig_dir, paste0(disease_id, "_single_filter_JI_heatmap.png")),
    p_ji,
    width = 210,
    height = 180,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

# ==============================================================================
# 5. Stability and conditioning metrics for single filters
# ==============================================================================

message("\n[metrics] Computing eRank, VIF and complementary stability metrics.")

single_metrics <- purrr::map_dfr(train_paths, function(task_path) {
  task_id <- clean_task_id(task_path)
  hc_label <- clean_hc_label(task_id, disease_id)

  df <- load_task_dataframe(task_path)

  task_feature_dir <- file.path(feature_dir, task_id)
  sets <- load_feature_sets_for_task(task_feature_dir)

  if (is.null(sets) || length(sets) == 0) return(NULL)

  feature_space_metrics(sets, df, n_boot = n_boot) |>
    mutate(
      task_id = task_id,
      HC = hc_label
    )
})

if (nrow(single_metrics) > 0) {
  single_metrics <- single_metrics |>
    filter(Filter %in% filter_ids) |>
    mutate(
      Method = dplyr::recode(as.character(Filter), !!!name_map),
      HC = factor(HC, levels = hc_order),
      Method = factor(Method, levels = unname(name_map))
    ) |>
    select(
      task_id,
      HC,
      Method,
      Filter,
      n_features,
      n_used,
      mean_vif,
      max_vif,
      MeanAbsCorr,
      Jaccard,
      eRank,
      eRank_lo,
      eRank_hi
    ) |>
    arrange(HC, Method)

  write_csv(
    single_metrics,
    file.path(table_dir, paste0(disease_id, "_single_filter_stability_conditioning.csv"))
  )
}

# ==============================================================================
# 6. Optional HybridA / HybridB stability and conditioning
# ==============================================================================

if (dir.exists(hybrid_dir)) {
  message("\n[hybrids] Hybrid directory found. Computing HybridA/HybridB metrics.")

  hybrid_metrics <- purrr::map_dfr(train_paths, function(task_path) {
    task_id <- clean_task_id(task_path)
    hc_label <- clean_hc_label(task_id, disease_id)

    df <- load_task_dataframe(task_path)

    fA <- file.path(hybrid_dir, paste0("hybrid_features_HybridA_", task_id, ".rds"))
    fB <- file.path(hybrid_dir, paste0("hybrid_features_HybridB_", task_id, ".rds"))

    sets <- list(
      HybridA = safe_read_vec(fA),
      HybridB = safe_read_vec(fB)
    )

    sets <- sets[lengths(sets) >= 2]

    if (length(sets) == 0) return(NULL)

    feature_space_metrics(sets, df, n_boot = n_boot) |>
      mutate(
        task_id = task_id,
        HC = hc_label,
        Method = Filter
      )
  })

  if (nrow(hybrid_metrics) > 0) {
    hybrid_metrics <- hybrid_metrics |>
      mutate(
        HC = factor(HC, levels = hc_order),
        Method = factor(Method, levels = c("HybridA", "HybridB"))
      ) |>
      select(
        task_id,
        HC,
        Method,
        Filter,
        n_features,
        n_used,
        mean_vif,
        max_vif,
        MeanAbsCorr,
        Jaccard,
        eRank,
        eRank_lo,
        eRank_hi
      ) |>
      arrange(HC, Method)

    write_csv(
      hybrid_metrics,
      file.path(table_dir, paste0(disease_id, "_hybrid_filter_stability_conditioning.csv"))
    )
  } else {
    message("[hybrids] No usable HybridA/HybridB files found.")
  }
} else {
  message("\n[hybrids] No hybrid directory found. Skipping HybridA/HybridB metrics.")
}

# ==============================================================================
# 7. Save reproducibility log
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(table_dir, paste0(disease_id, "_sessionInfo.txt"))
)

message("\n[done] Script 02 completed successfully.")
message("[done] Feature subsets written to: ", feature_dir)
message("[done] Tables written to: ", table_dir)
