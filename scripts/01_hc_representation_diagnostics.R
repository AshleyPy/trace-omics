# ==============================================================================
# TRACE-Omics
# Script 01: Harmonisation-configuration comparisons and representation diagnostics
#
# Purpose:
#   This script compares RNA-seq harmonisation configurations and evaluates
#   representation quality before feature selection and machine-learning modelling.
#
# Required input structure:
#
#   data/raw/TB/
#   ├── counts/
#   │   ├── cohort_1.csv
#   │   ├── cohort_2.csv
#   │   └── cohort_3.csv
#   └── metadata.rds
#
# Count-table requirement:
#   - Each count table must contain a column called: symbol
#   - Remaining columns must be sample IDs
#
# Metadata requirement:
#   - Must contain: sample_id, cohort, outcome
#
# Main outputs:
#   - Harmonised matrices
#   - QC tables
#   - Depth-dependence summaries
#   - PCA plots
#   - Silhouette summaries
#   - sessionInfo log
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(DESeq2)
  library(sva)
  library(limma)
  library(ggplot2)
  library(scales)
  library(cluster)
  library(RColorBrewer)
})

set.seed(7832)

# ==============================================================================
# 1. User settings
# ==============================================================================

disease_id <- "TB"

raw_dir     <- file.path("data", "raw", disease_id)
counts_dir  <- file.path(raw_dir, "counts")
metadata_fp <- file.path(raw_dir, "metadata.rds")

out_dir <- file.path("results", "01_hc_representation_diagnostics", disease_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_col  <- "sample_id"
cohort_col  <- "cohort"
outcome_col <- "outcome"

min_gene_count <- 10
reference_batch <- NULL

method_levels <- c(
  "RAW_counts",
  "VST_only",
  "VST_ComBat",
  "VST_ComBatRef",
  "VST_SVA",
  "ComBatSeq_log1p"
)

method_labels <- c(
  RAW_counts      = "RAW counts",
  VST_only        = "VST only",
  VST_ComBat      = "VST + ComBat",
  VST_ComBatRef   = "VST + ComBatRef",
  VST_SVA         = "VST + SVA",
  ComBatSeq_log1p = "ComBat-Seq + log1p"
)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

read_counts <- function(fp) {
  message("[read] ", fp)

  df <- suppressMessages(read_csv(fp, show_col_types = FALSE))

  if (!"symbol" %in% colnames(df)) {
    stop("Count table must contain a 'symbol' column: ", fp, call. = FALSE)
  }

  df <- df |>
    filter(!is.na(symbol), nzchar(symbol)) |>
    mutate(symbol = str_trim(symbol))

  df <- df[!duplicated(df$symbol), ]

  df
}

merge_count_tables <- function(counts_dir) {
  count_files <- list.files(counts_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(count_files) < 2) {
    stop("At least two cohort count tables are required in: ", counts_dir, call. = FALSE)
  }

  counts_list <- lapply(count_files, read_counts)
  names(counts_list) <- tools::file_path_sans_ext(basename(count_files))

  shared_genes <- Reduce(intersect, lapply(counts_list, \(x) x$symbol))

  if (length(shared_genes) < 1000) {
    warning("Fewer than 1000 shared genes detected across cohorts.")
  }

  message("[merge] shared genes across cohorts: ", length(shared_genes))

  counts_list <- lapply(counts_list, function(df) {
    df |>
      filter(symbol %in% shared_genes) |>
      arrange(symbol)
  })

  merged <- Reduce(
    function(x, y) inner_join(x, y, by = "symbol"),
    counts_list
  )

  merged <- merged |>
    mutate(symbol = make.unique(symbol)) |>
    column_to_rownames("symbol")

  mat <- as.matrix(merged)
  storage.mode(mat) <- "numeric"

  if (anyNA(mat)) stop("Merged count matrix contains NA values.", call. = FALSE)
  if (!all(is.finite(mat))) stop("Merged count matrix contains non-finite values.", call. = FALSE)

  mat <- round(mat)
  storage.mode(mat) <- "integer"

  mat
}

read_metadata <- function(metadata_fp, sample_col, cohort_col, outcome_col) {
  if (!file.exists(metadata_fp)) {
    stop("Metadata file not found: ", metadata_fp, call. = FALSE)
  }

  meta <- as.data.table(readRDS(metadata_fp))

  required <- c(sample_col, cohort_col, outcome_col)
  missing <- setdiff(required, names(meta))

  if (length(missing) > 0) {
    stop("Metadata is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  meta[, (sample_col) := trimws(as.character(get(sample_col)))]
  meta[, (cohort_col) := droplevels(factor(get(cohort_col)))]
  meta[, (outcome_col) := droplevels(factor(get(outcome_col)))]

  if (anyDuplicated(meta[[sample_col]]) > 0) {
    stop("Metadata contains duplicated sample IDs.", call. = FALSE)
  }

  meta
}

align_counts_metadata <- function(counts, meta, sample_col) {
  colnames(counts) <- trimws(as.character(colnames(counts)))

  shared_samples <- intersect(colnames(counts), meta[[sample_col]])

  if (length(shared_samples) < 10) {
    stop("Very few samples matched between counts and metadata. Check sample IDs.",
         call. = FALSE)
  }

  message("[align] matched samples: ", length(shared_samples))
  message("[align] count samples dropped: ", length(setdiff(colnames(counts), shared_samples)))
  message("[align] metadata samples dropped: ", length(setdiff(meta[[sample_col]], shared_samples)))

  counts <- counts[, shared_samples, drop = FALSE]
  meta <- meta[match(shared_samples, meta[[sample_col]]), ]

  stopifnot(identical(colnames(counts), meta[[sample_col]]))

  list(counts = counts, meta = meta)
}

coerce_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  if (anyNA(x)) stop("Matrix contains NA values.", call. = FALSE)
  if (!all(is.finite(x))) stop("Matrix contains non-finite values.", call. = FALSE)

  x
}

align_matrix_list <- function(mat_list, meta, sample_col) {
  mat_list <- lapply(mat_list, coerce_numeric_matrix)

  shared_samples <- Reduce(intersect, c(lapply(mat_list, colnames), list(meta[[sample_col]])))

  if (length(shared_samples) < 10) {
    stop("Too few shared samples across harmonised matrices.", call. = FALSE)
  }

  mat_list <- lapply(mat_list, function(m) m[, shared_samples, drop = FALSE])
  meta <- meta[match(shared_samples, meta[[sample_col]]), ]

  for (nm in names(mat_list)) {
    stopifnot(identical(colnames(mat_list[[nm]]), meta[[sample_col]]))
  }

  list(mat_list = mat_list, meta = meta)
}

pick_detect_threshold <- function(mat) {
  vals <- mat[is.finite(mat)]
  integerish <- all(abs(vals - round(vals)) < 1e-8)

  if (integerish && min(vals, na.rm = TRUE) >= 0) {
    return(1)
  }

  1e-8
}

qc_from_matrix <- function(mat, method_name, meta, sample_col, cohort_col, outcome_col) {
  thr <- pick_detect_threshold(mat)

  dt <- data.table(
    sample = colnames(mat),
    total_signal = colSums(mat, na.rm = TRUE),
    nonzero_features = colSums(mat > thr, na.rm = TRUE),
    method = method_name
  )

  meta_small <- meta[, .(
    sample = get(sample_col),
    cohort = get(cohort_col),
    outcome = get(outcome_col)
  )]

  merge(dt, meta_small, by = "sample", all.x = TRUE, sort = FALSE)
}

sanitize_for_pca <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  mat <- mat[rowSums(!is.finite(mat)) == 0, , drop = FALSE]

  gene_var <- apply(mat, 1, var)
  mat <- mat[is.finite(gene_var) & gene_var > 0, , drop = FALSE]

  mat
}

pca_scores <- function(mat, method_name) {
  mat <- sanitize_for_pca(mat)

  pr <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_explained <- (pr$sdev^2) / sum(pr$sdev^2) * 100

  data.table(
    sample = rownames(pr$x),
    PC1 = pr$x[, 1],
    PC2 = pr$x[, 2],
    method = method_name,
    PC1_var = var_explained[1],
    PC2_var = var_explained[2]
  )
}

make_hulls <- function(scores_dt) {
  scores_dt[, {
    if (.N >= 3) .SD[grDevices::chull(PC1, PC2)] else .SD[0]
  }, by = .(method, cohort)]
}

run_pca_full <- function(mat) {
  mat <- sanitize_for_pca(mat)
  pr <- prcomp(t(mat), center = TRUE, scale. = TRUE)

  list(
    scores = pr$x,
    var = (pr$sdev^2) / sum(pr$sdev^2)
  )
}

silhouette_over_k <- function(scores, labels, k_grid, var_vec, var_cap = 0.90) {
  labels <- droplevels(factor(labels))

  if (length(levels(labels)) < 2) {
    return(data.table())
  }

  cumvar <- cumsum(var_vec)

  k_grid <- k_grid[k_grid <= ncol(scores)]
  k_grid <- k_grid[cumvar[k_grid] <= var_cap | k_grid == min(k_grid)]

  if (length(k_grid) == 0) {
    k_grid <- min(2, ncol(scores))
  }

  rbindlist(lapply(k_grid, function(k) {
    xk <- scores[, seq_len(k), drop = FALSE]
    d <- dist(xk)

    sil <- silhouette(as.integer(labels), d)

    data.table(
      k = k,
      mean_sil = mean(sil[, "sil_width"]),
      median_sil = median(sil[, "sil_width"]),
      cumvar = cumvar[k]
    )
  }))
}

choose_k_star <- function(dt) {
  max_s <- max(dt$mean_sil, na.rm = TRUE)
  dt[mean_sil >= 0.95 * max_s][order(k)][1]
}

theme_trace <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.border = element_rect(colour = "black", linewidth = 0.4)
  )

# ==============================================================================
# 3. Load, merge, and align data
# ==============================================================================

counts <- merge_count_tables(counts_dir)
meta <- read_metadata(metadata_fp, sample_col, cohort_col, outcome_col)

aligned <- align_counts_metadata(counts, meta, sample_col)
counts <- aligned$counts
meta <- aligned$meta

message("[input] final count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

saveRDS(counts, file.path(out_dir, paste0(disease_id, "_joined_counts_matrix.rds")))
saveRDS(meta, file.path(out_dir, paste0(disease_id, "_aligned_metadata.rds")))

# ==============================================================================
# 4. Harmonisation configurations
# ==============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = as.data.frame(meta),
  design = as.formula(paste("~", cohort_col, "+", outcome_col))
)

keep <- rowSums(counts(dds)) >= min_gene_count
dds <- dds[keep, ]

vsd <- suppressMessages(varianceStabilizingTransformation(dds, blind = TRUE))
expr_vst <- assay(vsd)

counts_filtered <- counts[rownames(expr_vst), , drop = FALSE]

mod_outcome <- model.matrix(
  as.formula(paste("~", outcome_col)),
  data = as.data.frame(meta)
)

vst_combat <- ComBat(
  dat = as.matrix(expr_vst),
  batch = meta[[cohort_col]],
  mod = mod_outcome,
  par.prior = TRUE,
  prior.plots = FALSE
)

if (is.null(reference_batch)) {
  reference_batch <- names(sort(table(meta[[cohort_col]]), decreasing = TRUE))[1]
}

vst_combat_ref <- ComBat(
  dat = as.matrix(expr_vst),
  batch = meta[[cohort_col]],
  mod = mod_outcome,
  par.prior = TRUE,
  prior.plots = FALSE,
  ref.batch = reference_batch
)

mod <- model.matrix(as.formula(paste("~", outcome_col)), data = as.data.frame(meta))
mod0 <- model.matrix(~ 1, data = as.data.frame(meta))

n_sv_be <- num.sv(as.matrix(expr_vst), mod, method = "be")
n_sv_leek <- num.sv(as.matrix(expr_vst), mod, method = "leek")

max_sv_allowed <- max(0, ncol(expr_vst) - ncol(mod) - 1)
n_sv_used <- max(0, min(n_sv_be, floor(ncol(expr_vst) / 3), max_sv_allowed))

if (n_sv_used > 0) {
  svobj <- sva(as.matrix(expr_vst), mod = mod, mod0 = mod0, n.sv = n_sv_used)

  expr_vst_sva <- removeBatchEffect(
    x = as.matrix(expr_vst),
    covariates = svobj$sv,
    design = mod
  )
} else {
  warning("SVA estimated 0 usable surrogate variables. VST_SVA is set equal to VST_only.")
  expr_vst_sva <- expr_vst
}

counts_cbseq <- counts_filtered
storage.mode(counts_cbseq) <- "integer"

cbseq_counts <- ComBat_seq(
  counts = counts_cbseq,
  batch = meta[[cohort_col]],
  group = meta[[outcome_col]]
)

cbseq_log1p <- log1p(cbseq_counts)

mat_list <- list(
  RAW_counts = counts_filtered,
  VST_only = expr_vst,
  VST_ComBat = vst_combat,
  VST_ComBatRef = vst_combat_ref,
  VST_SVA = expr_vst_sva,
  ComBatSeq_log1p = cbseq_log1p
)

aligned_hc <- align_matrix_list(mat_list, meta, sample_col)
mat_list <- aligned_hc$mat_list
meta <- aligned_hc$meta

saveRDS(mat_list, file.path(out_dir, paste0(disease_id, "_harmonised_matrices.rds")))

saveRDS(
  list(
    reference_batch = reference_batch,
    n_sv_be = n_sv_be,
    n_sv_leek = n_sv_leek,
    n_sv_used = n_sv_used,
    retained_genes = rownames(expr_vst)
  ),
  file.path(out_dir, paste0(disease_id, "_harmonisation_parameters.rds"))
)

# ==============================================================================
# 5. QC diagnostics
# ==============================================================================

qc_all <- rbindlist(
  Map(
    qc_from_matrix,
    mat = mat_list,
    method_name = names(mat_list),
    MoreArgs = list(
      meta = meta,
      sample_col = sample_col,
      cohort_col = cohort_col,
      outcome_col = outcome_col
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

raw_counts <- as.matrix(mat_list$RAW_counts)

true_detected_raw <- data.table(
  sample = colnames(raw_counts),
  detected_genes_raw = colSums(raw_counts > 0, na.rm = TRUE)
)

qc_all <- merge(qc_all, true_detected_raw, by = "sample", all.x = TRUE)

qc_all[, method := factor(method, levels = method_levels)]
qc_all[, cohort := droplevels(factor(cohort))]
qc_all[, outcome := droplevels(factor(outcome))]

summary_qc <- qc_all[, .(
  n_samples = .N,
  median_total = median(total_signal, na.rm = TRUE),
  IQR_total = IQR(total_signal, na.rm = TRUE),
  median_nonzero = median(nonzero_features, na.rm = TRUE),
  IQR_nonzero = IQR(nonzero_features, na.rm = TRUE),
  median_detected_raw = median(detected_genes_raw, na.rm = TRUE),
  IQR_detected_raw = IQR(detected_genes_raw, na.rm = TRUE)
), by = .(cohort, method)][order(cohort, method)]

fwrite(summary_qc, file.path(out_dir, paste0(disease_id, "_qc_summary.csv")))
saveRDS(qc_all, file.path(out_dir, paste0(disease_id, "_qc_sample_level.rds")))

# ==============================================================================
# 6. Depth-dependence diagnostics
# ==============================================================================

corr_overall <- qc_all[, .(
  rho = cor(
    detected_genes_raw,
    total_signal,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
), by = method][order(method)]

corr_by_cohort <- qc_all[, .(
  rho = cor(
    detected_genes_raw,
    total_signal,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
), by = .(method, cohort)][order(method, cohort)]

fwrite(corr_overall, file.path(out_dir, paste0(disease_id, "_depth_dependence_overall.csv")))
fwrite(corr_by_cohort, file.path(out_dir, paste0(disease_id, "_depth_dependence_by_cohort.csv")))

# ==============================================================================
# 7. QC plots
# ==============================================================================

geom_violin_box <- list(
  geom_violin(trim = FALSE, linewidth = 0.3, colour = "black"),
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.4, colour = "black"),
  stat_summary(fun = median, geom = "point", size = 1.2)
)

p_total <- ggplot(qc_all, aes(x = cohort, y = total_signal)) +
  geom_violin_box +
  facet_wrap(~ method, scales = "free_y", labeller = as_labeller(method_labels)) +
  scale_y_continuous(labels = label_number(scale_cut = cut_si("unit"))) +
  labs(
    title = paste0(disease_id, ": total signal across harmonisation configurations"),
    x = "Cohort",
    y = "Total signal"
  ) +
  theme_trace

p_nonzero <- ggplot(qc_all, aes(x = cohort, y = nonzero_features)) +
  geom_violin_box +
  facet_wrap(~ method, scales = "free_y", labeller = as_labeller(method_labels)) +
  labs(
    title = paste0(disease_id, ": non-zero features across harmonisation configurations"),
    x = "Cohort",
    y = "Non-zero features"
  ) +
  theme_trace

p_detected_raw <- ggplot(qc_all[method == "RAW_counts"], aes(x = cohort, y = detected_genes_raw)) +
  geom_violin_box +
  labs(
    title = paste0(disease_id, ": detected genes in raw counts"),
    x = "Cohort",
    y = "Detected genes"
  ) +
  theme_trace

ggsave(file.path(out_dir, paste0(disease_id, "_qc_total_signal.png")),
       p_total, width = 180, height = 120, units = "mm", dpi = 300)

ggsave(file.path(out_dir, paste0(disease_id, "_qc_nonzero_features.png")),
       p_nonzero, width = 180, height = 120, units = "mm", dpi = 300)

ggsave(file.path(out_dir, paste0(disease_id, "_qc_detected_genes_raw.png")),
       p_detected_raw, width = 120, height = 100, units = "mm", dpi = 300)

# ==============================================================================
# 8. PCA diagnostics
# ==============================================================================

pca_dt_list <- Map(pca_scores, mat_list, names(mat_list))
scores_all <- rbindlist(pca_dt_list, use.names = TRUE, fill = TRUE)

meta_small <- meta[, .(
  sample = get(sample_col),
  cohort = get(cohort_col),
  outcome = get(outcome_col)
)]

scores_all <- merge(scores_all, meta_small, by = "sample", all.x = TRUE)
scores_all[, method := factor(method, levels = method_levels)]

var_tbl <- scores_all[, .(
  PC1_var = unique(PC1_var)[1],
  PC2_var = unique(PC2_var)[1]
), by = method][order(method)]

fwrite(scores_all, file.path(out_dir, paste0(disease_id, "_pca_scores.csv")))
fwrite(var_tbl, file.path(out_dir, paste0(disease_id, "_pca_variance.csv")))

hulls_all <- make_hulls(scores_all)

xlim_all <- range(scores_all$PC1, na.rm = TRUE)
ylim_all <- range(scores_all$PC2, na.rm = TRUE)

pad_range <- function(x, p = 0.05) {
  width <- diff(x)
  c(x[1] - p * width, x[2] + p * width)
}

xlim_all <- pad_range(xlim_all)
ylim_all <- pad_range(ylim_all)

strip_labs <- setNames(
  paste0(
    method_labels[method_levels],
    "\nPC1 ",
    sprintf("%.1f", var_tbl$PC1_var[match(method_levels, var_tbl$method)]),
    "%, PC2 ",
    sprintf("%.1f", var_tbl$PC2_var[match(method_levels, var_tbl$method)]),
    "%"
  ),
  method_levels
)

cohort_levels <- levels(droplevels(factor(scores_all$cohort)))
palette_set1 <- brewer.pal(max(3, length(cohort_levels)), "Set1")
cohort_cols <- setNames(palette_set1[seq_along(cohort_levels)], cohort_levels)

p_pca <- ggplot() +
  geom_polygon(
    data = hulls_all,
    aes(PC1, PC2, group = interaction(method, cohort), fill = cohort),
    alpha = 0.15,
    colour = NA
  ) +
  geom_polygon(
    data = hulls_all,
    aes(PC1, PC2, group = interaction(method, cohort), colour = cohort),
    fill = NA,
    linewidth = 0.4
  ) +
  geom_point(
    data = scores_all,
    aes(PC1, PC2, colour = cohort, shape = outcome),
    size = 1.7,
    alpha = 0.9
  ) +
  facet_wrap(~ method, ncol = 3, labeller = as_labeller(strip_labs)) +
  scale_colour_manual(values = cohort_cols) +
  scale_fill_manual(values = cohort_cols) +
  coord_fixed() +
  scale_x_continuous(limits = xlim_all) +
  scale_y_continuous(limits = ylim_all) +
  labs(
    title = paste0(disease_id, ": comparable PCA across harmonisation configurations"),
    x = "PC1",
    y = "PC2",
    colour = "Cohort",
    fill = "Cohort",
    shape = "Outcome"
  ) +
  theme_trace +
  theme(panel.grid = element_blank())

ggsave(file.path(out_dir, paste0(disease_id, "_pca_hc_comparison.png")),
       p_pca, width = 180, height = 130, units = "mm", dpi = 300)

# ==============================================================================
# 9. Silhouette diagnostics
# ==============================================================================

pca_objs <- lapply(mat_list, run_pca_full)
k_grid <- 2:50

sil_results <- rbindlist(lapply(names(pca_objs), function(meth) {
  scores <- pca_objs[[meth]]$scores
  varvec <- pca_objs[[meth]]$var

  scores <- scores[meta[[sample_col]], , drop = FALSE]

  sil_cohort <- silhouette_over_k(scores, meta[[cohort_col]], k_grid, varvec)
  if (nrow(sil_cohort) > 0) {
    sil_cohort[, `:=`(method = meth, label = "cohort")]
  }

  sil_outcome <- silhouette_over_k(scores, meta[[outcome_col]], k_grid, varvec)
  if (nrow(sil_outcome) > 0) {
    sil_outcome[, `:=`(method = meth, label = "outcome")]
  }

  rbind(sil_cohort, sil_outcome, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

sil_results[, method := factor(method, levels = method_levels)]

k_star <- sil_results[, choose_k_star(.SD), by = .(method, label)]

setnames(
  k_star,
  old = c("k", "mean_sil", "median_sil", "cumvar"),
  new = c("k_star", "mean_sil_at_k", "median_sil_at_k", "cumvar_at_k")
)

sil_summary <- k_star[, .(
  method,
  label,
  k_star,
  mean_sil = mean_sil_at_k,
  median_sil = median_sil_at_k,
  cumvar = cumvar_at_k
)][order(method, label)]

fwrite(sil_results, file.path(out_dir, paste0(disease_id, "_silhouette_all_k.csv")))
fwrite(sil_summary, file.path(out_dir, paste0(disease_id, "_silhouette_summary.csv")))

sil_cols <- c(
  cohort = brewer.pal(3, "Set1")[1],
  outcome = brewer.pal(3, "Set1")[2]
)

p_sil <- ggplot(sil_results, aes(x = k, y = mean_sil, colour = label)) +
  geom_line(linewidth = 0.7, alpha = 0.95) +
  geom_point(size = 1.1, alpha = 0.9) +
  geom_vline(
    data = k_star,
    aes(xintercept = k_star, colour = label),
    linetype = "22",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  facet_wrap(~ method, ncol = 3, labeller = as_labeller(method_labels)) +
  scale_colour_manual(values = sil_cols, name = "Label") +
  coord_cartesian(ylim = c(-0.25, 0.95)) +
  labs(
    title = paste0(disease_id, ": silhouette profiles across principal components"),
    x = "Number of PCs",
    y = "Mean silhouette"
  ) +
  theme_trace +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, paste0(disease_id, "_silhouette_profiles.png")),
       p_sil, width = 180, height = 120, units = "mm", dpi = 300)

# ==============================================================================
# 10. Save reproducibility log
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, paste0(disease_id, "_sessionInfo.txt"))
)

message("[done] Script 01 completed successfully.")
message("[done] Outputs written to: ", out_dir)
