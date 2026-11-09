# TRACE-Omics | Core configuration, validation and orchestration
# Designed for binary, gene-level bulk RNA-seq tasks from any disease.
# No patient-level data are distributed with this software.
# This refactor is a publication candidate, not a validated numerical reproduction.

trace_seed <- function(..., base = 7832L) {
  keys <- as.numeric(unlist(list(...), use.names = FALSE))
  if (any(!is.finite(keys))) stop("Subseed keys must be finite numbers.")
  z <- (as.double(base) + sum(keys * (1009 * seq_along(keys)))) %% 2000000000
  as.integer(z)
}

trace_rng <- function(seed = 7832L) {
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max) stop("Invalid non-negative integer RNG seed.")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  invisible(seed)
}

trace_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required package(s): ", paste(absent, collapse = ", "))
  invisible(TRUE)
}

trace_find_root <- function(config_file = NULL, project_root = NULL) {
  if (!is.null(project_root)) {
    root <- normalizePath(path.expand(project_root), winslash = "/", mustWork = TRUE)
    if (!file.exists(file.path(root, "R", "supplied_functions.R")))
      stop("Specified project_root does not contain R/supplied_functions.R: ", root)
    return(root)
  }
  starts <- c(getwd(), if (!is.null(config_file)) dirname(normalizePath(config_file)) else character())
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(current, "R", "supplied_functions.R"))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate TRACE-Omics project root. Supply project_root explicitly.")
}

trace_path <- function(path, root) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid path setting.")
  path <- path.expand(path)
  # Covers Windows drive paths, UNC paths and Unix absolute paths.
  if (grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2}|/)", path)) path else file.path(root, path)
}

trace_read_config <- function(config_file) {
  if (!file.exists(config_file)) stop("Configuration not found: ", config_file)
  extension <- tolower(tools::file_ext(config_file))
  if (extension %in% c("yml", "yaml")) {
    trace_require("yaml")
    cfg <- yaml::read_yaml(config_file)
  } else if (extension == "rds") {
    cfg <- readRDS(config_file)
  } else {
    stop("Configuration must be YAML or RDS. Do not use executable dget() files.")
  }
  if (!is.list(cfg)) stop("Configuration must contain a named list.")
  cfg
}

trace_read_counts <- function(path, gene_column = "gene_id") {
  if (!file.exists(path)) stop("Count file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
  } else if (ext %in% c("csv", "tsv")) {
    tab <- if (ext == "tsv") utils::read.delim(path, check.names = FALSE) else
      utils::read.csv(path, check.names = FALSE)
    if (!gene_column %in% names(tab)) stop("Count table requires gene column: ", gene_column)
    genes <- as.character(tab[[gene_column]])
    x <- as.matrix(tab[, setdiff(names(tab), gene_column), drop = FALSE])
    rownames(x) <- genes
  } else {
    stop("Counts must be .rds, .csv or .tsv.")
  }
  if (!is.matrix(x) || is.null(rownames(x)) || is.null(colnames(x)))
    stop("Counts must be a named genes-by-samples matrix.")
  if (!is.numeric(x)) {
    original <- x
    suppressWarnings(storage.mode(x) <- "numeric")
    if (anyNA(x) && !anyNA(original)) stop("Count table contains non-numeric values.")
  }
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0) ||
      any(abs(x - round(x)) > 1e-8)) {
    stop("Expected unnormalised, finite, non-negative integer gene counts.")
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x)))) {
    stop("Gene and sample identifiers must be non-empty and unique.")
  }
  if (any(colSums(x) == 0)) stop("At least one sample has zero total counts.")
  x
}

trace_read_metadata <- function(path, columns = list(),
                                assume_one_sample_per_participant = FALSE) {
  if (!file.exists(path)) stop("Metadata not found: ", path)
  ext <- tolower(tools::file_ext(path))
  m <- if (ext == "rds") readRDS(path) else if (ext == "csv")
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE) else if (ext == "tsv")
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE) else
      stop("Metadata must be .rds, .csv or .tsv.")
  m <- as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "study", "outcome")
  source_columns <- vapply(required, function(k) {
    z <- columns[[k]]
    if (is.null(z)) k else as.character(z)
  }, character(1))
  absent <- setdiff(source_columns, names(m))
  if (length(absent)) stop("Metadata missing columns: ", paste(absent, collapse = ", "))
  for (k in required) m[[k]] <- trimws(as.character(m[[source_columns[[k]]]]))
  participant_source <- if (is.null(columns$participant_id)) "participant_id" else columns$participant_id
  if (participant_source %in% names(m)) {
    m$participant_id <- trimws(as.character(m[[participant_source]]))
  } else if (isTRUE(assume_one_sample_per_participant)) {
    warning("Using sample_id as participant_id. Verify there is only one sample per participant.")
    m$participant_id <- m$sample_id
  } else {
    stop("Metadata require participant_id unless one_sample_per_participant is explicitly TRUE.")
  }
  if (anyNA(m[, c(required, "participant_id")]) ||
      any(!nzchar(as.matrix(m[, c(required, "participant_id")])))) {
    stop("Required metadata fields contain missing or blank values.")
  }
  if (anyDuplicated(m$sample_id)) stop("Duplicated sample identifiers in metadata.")
  m
}

trace_make_partition <- function(metadata, fraction = 0.8, seed = 7832L,
                                 strata = c("study", "outcome")) {
  if (length(fraction) != 1L || fraction <= 0 || fraction >= 1)
    stop("development_fraction must be strictly between 0 and 1.")
  trace_rng(seed)
  # Assign participants, not samples, to prevent repeated-measure leakage.
  id_study <- split(metadata$study, metadata$participant_id)
  id_outcome <- split(metadata$outcome, metadata$participant_id)
  if (any(lengths(lapply(id_study, unique)) != 1L) ||
      any(lengths(lapply(id_outcome, unique)) != 1L))
    stop("A participant spans study or outcome strata; provide an explicit validated manifest.")
  u <- metadata[!duplicated(metadata$participant_id), c("participant_id", "study", "outcome")]
  if (length(strata) < 1L || !all(strata %in% names(u)))
    stop("Invalid split strata. Choose from the columns study and outcome.")
  group <- do.call(interaction, c(unname(u[, strata, drop = FALSE]), list(drop = TRUE)))
  partition <- rep("holdout", nrow(u))
  for (ids in split(seq_len(nrow(u)), group)) {
    if (length(ids) < 2L) stop("Cannot split a study-outcome stratum with fewer than two participants.")
    ndev <- min(length(ids) - 1L, max(1L, round(fraction * length(ids))))
    partition[sample(ids, ndev)] <- "development"
  }
  result <- data.frame(sample_id = metadata$sample_id,
                       partition = partition[match(metadata$participant_id, u$participant_id)])
  stopifnot(!anyNA(result$partition))
  result
}

trace_load <- function(config_file) {
  cfg <- trace_read_config(config_file)
  required <- c("counts", "metadata", "output", "seed")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Configuration missing: ", paste(missing, collapse = ", "))
  root <- trace_find_root(config_file, cfg$project_root)
  inputs <- lapply(cfg[c("counts", "metadata")], trace_path, root = root)
  x <- trace_read_counts(inputs$counts, gene_column = if (is.null(cfg$gene_column)) "gene_id" else cfg$gene_column)
  m <- trace_read_metadata(inputs$metadata,
      columns = if (is.null(cfg$columns)) list() else cfg$columns,
      assume_one_sample_per_participant = isTRUE(cfg$one_sample_per_participant))
  class_map <- cfg$class_map
  if (!is.null(class_map)) {
    if (!all(c("case", "control") %in% names(class_map)))
      stop("class_map requires case and control entries.")
    positives <- as.character(class_map$case)
    negatives <- as.character(class_map$control)
    if (length(intersect(positives, negatives))) stop("Class label mappings overlap.")
    m$outcome[m$outcome %in% positives] <- "case"
    m$outcome[m$outcome %in% negatives] <- "control"
  }
  if (!all(m$outcome %in% c("case", "control")))
    stop("Outcome must be case/control, or map original labels using class_map.")
  if (!setequal(colnames(x), m$sample_id))
    stop("Count matrix and metadata must contain identical sample sets. Review input harmonisation.")
  m <- m[match(colnames(x), m$sample_id), , drop = FALSE]
  if (!identical(colnames(x), m$sample_id)) stop("Sample alignment failure.")
  if (length(unique(m$outcome)) < 2L) stop("Both classes are required.")
  m$outcome <- factor(m$outcome, levels = c("control", "case"))
  m$disease <- m$outcome  # Compatibility with legacy functions (outcome, not disease name).
  m$cohort_group <- factor(m$study)
  m$cohort <- factor(m$study)
  mode <- if (is.null(cfg$partition_mode)) "manifest" else cfg$partition_mode
  if (mode == "manifest") {
    if (is.null(cfg$partition)) stop("partition manifest required for reproducibility mode.")
    partition_path <- trace_path(cfg$partition, root)
    z <- utils::read.csv(partition_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("sample_id", "partition") %in% names(z)) || anyDuplicated(z$sample_id) ||
        !setequal(z$sample_id, m$sample_id)) stop("Partition manifest must identify every sample exactly once.")
  } else if (mode == "generate") {
    z <- trace_make_partition(m,
        fraction = if (is.null(cfg$development_fraction)) 0.8 else cfg$development_fraction,
        seed = cfg$seed,
        strata = if (is.null(cfg$partition_strata)) c("study", "outcome") else cfg$partition_strata)
    partition_path <- NULL
  } else {
    stop("partition_mode must be manifest or generate.")
  }
  m$partition <- z$partition[match(m$sample_id, z$sample_id)]
  if (anyNA(m$partition) || !all(m$partition %in% c("development", "holdout")))
    stop("Partition labels must be development or holdout.")
  if (any(vapply(split(m$partition, m$participant_id), function(v) length(unique(v)) != 1L,
                 logical(1)))) stop("Participant crosses development/holdout partition.")
  for (part in c("development", "holdout"))
    if (length(unique(m$outcome[m$partition == part])) != 2L)
      stop("Both classes must appear in ", part, ".")
  # Keep experiment outputs separate from the repository itself.
  output <- trace_path(cfg$output, root)
  output <- normalizePath(output, winslash = "/", mustWork = FALSE)
  if (identical(output, root) || startsWith(output, paste0(root, "/")))
    stop("Output directory must be outside the Git repository to protect individual-level results.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  cfg$output <- output
  cfg$project_root <- root
  cfg$counts <- inputs$counts
  cfg$metadata <- inputs$metadata
  if (!is.null(partition_path)) cfg$partition <- partition_path
  trace_rng(cfg$seed)
  list(config = cfg, counts = x, metadata = m, partition_manifest = z)
}

trace_record <- function(d, stage) {
  if (!grepl("^[A-Za-z0-9_-]+$", stage)) stop("Stage name must be a simple identifier.")
  target <- file.path(d$config$output, stage)
  if (dir.exists(target)) stop("Stage output already exists: ", target)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(target, "sessionInfo.txt"))
  dput(d$config, file.path(target, "config.R"))
  saveRDS(list(kind = RNGkind(), seed = d$config$seed,
               state = get(".Random.seed", envir = .GlobalEnv)), file.path(target, "rng.rds"))
  input_paths <- unlist(d$config[intersect(c("counts", "metadata", "partition"), names(d$config))])
  input_paths <- input_paths[file.exists(input_paths)]
  utils::write.csv(data.frame(path = input_paths, md5 = unname(tools::md5sum(input_paths))),
                   file.path(target, "input_checksums.csv"), row.names = FALSE)
  utils::write.csv(d$partition_manifest, file.path(target, "partition_manifest.csv"), row.names = FALSE)
  invisible(target)
}

trace_source_functions <- function(project_root = NULL, seed = 7832L,
                                   min_total_training_count = 10L, top_n = 50L,
                                   max_sv = 50L) {
  root <- trace_find_root(project_root = project_root)
  e <- new.env(parent = globalenv())
  e$SEED <- as.integer(seed)
  e$POSITIVE_CLASS <- "case"
  e$NEGATIVE_CLASS <- "control"
  e$MIN_TOTAL_TRAINING_COUNT <- as.integer(min_total_training_count)
  e$TOP_N <- as.integer(top_n)
  e$MAX_SV <- as.integer(max_sv)
  sys.source(file.path(root, "R", "supplied_functions.R"), envir = e)
  e
}

trace_representation <- function(id, train, new = NULL, train_meta, new_meta = NULL,
                                 seed = 7832L, config = list(), project_root = NULL) {
  opts <- if (is.null(config$preprocessing)) list() else config$preprocessing
  get_opt <- function(name, default) if (is.null(opts[[name]])) default else opts[[name]]
  e <- trace_source_functions(project_root, seed,
      min_total_training_count = get_opt("minimum_total_count", 10L),
      top_n = get_opt("minimum_genes", 50L), max_sv = get_opt("max_surrogate_variables", 50L))
  if (id == "RAW_counts") return(list(train = train, new = new, state = list(genes = rownames(train))))
  if (id == "logCPM") {
    keep <- rowSums(train) >= get_opt("minimum_total_count", 10L)
    if (!any(keep)) stop("No genes survive training-only count filtering.")
    transform <- function(x) {
      lib <- colSums(x[keep, , drop = FALSE])
      if (any(lib <= 0)) stop("Zero library size after training gene filter.")
      log2(sweep(x[keep, , drop = FALSE], 2, lib, "/") * 1e6 + 1)
    }
    return(list(train = transform(train), new = if (!is.null(new)) transform(new) else NULL,
                state = list(genes = rownames(train)[keep])))
  }
  trace_require(c("DESeq2", "SummarizedExperiment", "matrixStats"))
  if (id == "VST_only") return(e$fit_apply_frozen_vst(train, new, train_meta))
  if (id == "VST_ComBat") {
    if (!is.null(new)) stop("Ordinary ComBat without a fitted new-data operator is descriptive only.")
    trace_require("neuroCombat")
    v <- e$fit_apply_frozen_vst(train, NULL, train_meta)
    fit <- neuroCombat::neuroCombat(dat = v$train, batch = as.character(train_meta$cohort_group),
                                     mod = NULL, parametric = TRUE)
    return(list(train = fit$dat.combat, new = NULL, state = fit$estimates))
  }
  if (id == "VST_fSVA") trace_require("sva")
  if (id %in% c("VST_ComBatRef", "ComBatSeq_log1p")) {
    if (!is.null(new) && is.null(new_meta)) stop(id, " requires new_meta with known cohorts.")
    if (!is.null(new_meta) && any(!new_meta$study %in% train_meta$study))
      stop(id, " cannot project an entirely unseen cohort using existing fitted batch parameters.")
    if (id == "VST_ComBatRef") trace_require("neuroCombat") else trace_require("edgeR")
  }
  e$fit_apply_representation(id, train, new, train_meta, new_meta, seed)
}

trace_select <- function(train, new, meta, method, k, seed, project_root = NULL) {
  if (length(k) != 1L || k < 1L || k > nrow(train)) stop("Invalid top-k feature count.")
  trace_require(c("mlr3", "mlr3filters", "mlr3learners", "matrixStats", "data.table"))
  e <- trace_source_functions(project_root, seed, top_n = k)
  e$select_top_features(train, new, meta, method, as.integer(k), seed)
}

trace_folds <- function(meta, k = 5L, seed = 7832L, manifest = NULL,
                        strata = c("study", "outcome")) {
  if (k < 2L) stop("Cross-validation needs at least two folds.")
  
    required_columns <- c(
      "sample_id",
      "participant_id",
      "study",
      "outcome"
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(meta)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "Missing required fold metadata: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (!is.null(manifest)) {
    z <- if (is.character(manifest)) utils::read.csv(manifest, stringsAsFactors = FALSE) else manifest
    if (!all(c("sample_id", "fold") %in% names(z)) || anyDuplicated(z$sample_id) ||
        !setequal(z$sample_id, meta$sample_id)) stop("Invalid frozen fold manifest.")
    f <- as.integer(z$fold[match(meta$sample_id, z$sample_id)])
    if (anyNA(f) || !setequal(unique(f), seq_len(k))) stop("Expected fold values 1..k.")
  } else {
    trace_rng(seed)
    u <- meta[!duplicated(meta$participant_id), , drop = FALSE]
    validate <- function(x) length(unique(x)) == 1L
    if (any(!vapply(split(meta$study, meta$participant_id), validate, logical(1))) ||
        any(!vapply(split(as.character(meta$outcome), meta$participant_id), validate, logical(1))))
      stop("Participants span strata; provide an explicit grouped fold manifest.")
    if (length(strata) < 1L || !all(strata %in% names(u)))
      stop("Invalid CV strata.")
    group <- do.call(interaction, c(unname(u[, strata, drop = FALSE]), list(drop = TRUE)))
    fg <- integer(nrow(u))
    for (idx in split(seq_len(nrow(u)), group)) {
      if (length(idx) < k) stop("Stratum has fewer participants than folds; reduce k or supply a manifest.")
      fg[sample(idx)] <- rep(seq_len(k), length.out = length(idx))
    }
    f <- fg[match(meta$participant_id, u$participant_id)]
  }
  if (any(vapply(split(f, meta$participant_id), function(z) length(unique(z)) > 1L, logical(1))))
    stop("Participant crosses CV folds.")
  for (j in seq_len(k)) {
    if (length(unique(meta$outcome[f == j])) != 2L ||
        length(unique(meta$outcome[f != j])) != 2L)
      stop("Missing class in CV train or validation fold ", j, ".")
  }
  f
}

trace_auc <- function(y, p) {
  trace_require("pROC")
  if (length(unique(y)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(as.character(y), p,
      levels = c("control", "case"), direction = "<", quiet = TRUE)))
}

trace_metrics <- function(y, p, threshold = 0.5) {
  trace_require(c("pROC", "PRROC"))
  if (length(y) != length(p) || anyNA(y) || anyNA(p) || any(!is.finite(p)) ||
      any(p < 0 | p > 1) || !all(y %in% c("case", "control")))
    stop("Invalid ground truth or probability predictions.")
  truth <- as.character(y) == "case"
  pred <- p >= threshold
  sens <- if (any(truth)) mean(pred[truth]) else NA_real_
  spec <- if (any(!truth)) mean(!pred[!truth]) else NA_real_
  tp <- sum(truth & pred); fp <- sum(!truth & pred); fn <- sum(truth & !pred)
  c(AUROC = trace_auc(y, p),
    AUPRC = if (length(unique(truth)) == 2L)
      PRROC::pr.curve(scores.class0 = p[truth], scores.class1 = p[!truth], curve = FALSE)$auc.integral else NA_real_,
    bACC = mean(c(sens, spec)), Sens = sens, Spec = spec,
    F1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_,
    Brier = mean((as.integer(truth) - p)^2), Prevalence = mean(truth))
}

trace_intervals <- function(y, p, threshold = 0.5, B = 2000L, seed = 7832L,
                            participant_id = NULL) {
  if (!is.null(participant_id) && anyDuplicated(participant_id))
    stop("Repeated participants require cluster-aware confidence intervals.")
  if (length(unique(y)) != 2L || B < 10L) stop("Both classes and B >= 10 are required.")
  point <- trace_metrics(y, p, threshold)
  point <- point[names(point) != "Prevalence"]
  trace_rng(seed)
  r <- pROC::roc(as.character(y), p, levels = c("control", "case"), direction = "<", quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  classes <- split(seq_along(y), as.character(y))
  boot <- replicate(B, {
    ii <- unlist(lapply(classes, function(z) sample(z, length(z), replace = TRUE)), use.names = FALSE)
    trace_metrics(y[ii], p[ii], threshold)[names(point)]
  })
  q <- t(apply(boot, 1L, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  q["AUROC", ] <- auc_ci[c(1L, 3L)]
  data.frame(metric = names(point), estimate = unname(point), lower = q[, 1L], upper = q[, 2L],
    method = ifelse(names(point) == "AUROC", "DeLong", "class-stratified percentile bootstrap"),
    row.names = NULL)
}

trace_calibration <- function(y, p, bins = 8L) {
  if (length(y) != length(p) || anyNA(y) || anyNA(p) || any(p < 0 | p > 1))
    stop("Invalid calibration inputs.")
  # Quantile-based bins; ties are combined, so the realised number may be smaller.
  breaks <- unique(stats::quantile(p, seq(0, 1, length.out = bins + 1L), names = FALSE))
  group <- if (length(breaks) < 2L) rep(1L, length(p)) else
    cut(p, breaks, include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(split(seq_along(p), group), function(i) {
    cases <- sum(y[i] == "case")
    ci <- stats::binom.test(cases, length(i))$conf.int
    data.frame(n = length(i), predicted = mean(p[i]), observed = cases / length(i),
               lower = ci[1L], upper = ci[2L])
  }))
}

trace_predict_pipeline <- function(train, new, meta, new_meta, spec, seed = 7832L,
                                   project_root = NULL, config = list()) {
  if (is.null(new) || is.null(new_meta)) stop("Independent evaluation data and metadata are required.")
  required <- c("representation", "filter", "learner", "genes")
  missing <- setdiff(required, names(spec))
  if (length(missing)) stop("Pipeline specification missing: ", paste(missing, collapse = ", "))
  e <- trace_source_functions(project_root, seed)
  r <- trace_representation(spec$representation, train, new, meta, new_meta, seed,
                             config = config, project_root = project_root)
  f <- trace_select(r$train, r$new, meta, spec$filter, as.integer(spec$genes), seed,
                     project_root = project_root)
  learner <- e$make_learner(spec$learner)
  if (!is.null(spec$parameters)) {
    for (name in names(spec$parameters)) learner$param_set$values[[name]] <- spec$parameters[[name]]
  }
  task <- e$make_classification_task(f$train, meta, "trace_training")
  trace_rng(seed)
  learner$train(task)
  x <- as.data.frame(t(f$validation), check.names = FALSE)
  names(x) <- make.names(rownames(f$validation), unique = TRUE)
  p <- as.numeric(learner$predict_newdata(x)$prob[, "case"])
  if (length(p) != nrow(new_meta) || anyNA(p) || any(p < 0 | p > 1))
    stop("Predicted probabilities failed validation.")
  list(probability = p, genes = f$genes, learner = learner,
       representation_state = r$state)
}
