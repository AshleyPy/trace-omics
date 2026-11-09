# Fast dependency-free contracts. These are NOT numerical reproduction tests.
source("R/core.R", local = FALSE)
source("R/utils_io.R", local = FALSE)
source("R/study_models.R", local = FALSE)
source("R/supplied_architecture_functions.R", local = FALSE)

stopifnot(identical(trace_lightgbm_grid()$Config_ID[1L], "LGBM_01"))
stopifnot(identical(clean_representation_label("VST + ComBat"), "VST + ComBat"))
stopifnot(identical(clean_representation_label("log2(CPM+1)"), "log2(CPM+1)"))
stopifnot(identical(clean_learner_label("C5.0"), "C5.0"))
stopifnot(identical(clean_feature_selection_label("auc"), "AUC"))

# 2 studies x 2 classes x 10 participants per stratum
m <- expand.grid(study = c("cohort_a", "cohort_b"), outcome = c("control", "case"),
                 sample_number = seq_len(10), stringsAsFactors = FALSE)
m$sample_id <- sprintf("sample_%02d", seq_len(nrow(m)))
m$participant_id <- m$sample_id
m$outcome <- factor(m$outcome, levels = c("control", "case"))
part <- trace_make_partition(m, fraction = 0.8, seed = 7832L)
stopifnot(nrow(part) == nrow(m), all(part$partition %in% c("development", "holdout")))
stopifnot(identical(part, trace_make_partition(m, fraction = 0.8, seed = 7832L)))
f <- trace_folds(m, k = 5L, seed = 7832L)
stopifnot(length(f) == nrow(m), setequal(unique(f), seq_len(5L)))
stopifnot(identical(f, trace_folds(m, k = 5L, seed = 7832L)))
cat("PASS base-R configuration, partition, fold and label checks.\n")

# End-to-end synthetic input/metadata/manifest loading, without ML packages.
tmp <- tempfile("trace-synthetic-")
dir.create(tmp)
set.seed(100)
counts <- matrix(as.integer(stats::rpois(200L * nrow(m), lambda = 20)),
                 nrow = 200L, ncol = nrow(m),
                 dimnames = list(sprintf("gene_%03d", 1:200), m$sample_id))
count_file <- file.path(tmp, "counts.rds")
metadata_file <- file.path(tmp, "metadata.csv")
partition_file <- file.path(tmp, "partition.csv")
config_file <- file.path(tmp, "config.rds")
saveRDS(counts, count_file)
utils::write.csv(m, metadata_file, row.names = FALSE)
utils::write.csv(part, partition_file, row.names = FALSE)
cfg <- list(counts = count_file, metadata = metadata_file, partition = partition_file,
            output = file.path(tmp, "private_output"), project_root = normalizePath(getwd()),
            seed = 7832L, partition_mode = "manifest")
saveRDS(cfg, config_file)
d <- trace_load(config_file)
stopifnot(identical(colnames(d$counts), d$metadata$sample_id),
          nrow(d$counts) == 200L,
          identical(sort(unique(d$metadata$partition)), c("development", "holdout")))
audit_dir <- trace_record(d, "synthetic_contract")
stopifnot(file.exists(file.path(audit_dir, "input_checksums.csv")),
          file.exists(file.path(audit_dir, "partition_manifest.csv")))
unlink(tmp, recursive = TRUE)
cat("PASS synthetic-data load and audit recording.\n")
