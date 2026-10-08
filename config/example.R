# Copy to config/local.R. Paths are relative to the repository root.
# Use an original frozen manifest for historical reproduction.
list(
 counts = "data/TB_counts.rds",
 metadata = "data/TB_metadata.csv",
 partition = "data/TB_partition_80_20.csv",
 output = "results/TB_run_001",
 seed = 7832L,
 representations = c("RAW_counts","VST_only","VST_fSVA","VST_ComBat","VST_ComBatRef","ComBatSeq_log1p"),
 filters = c("anova","auc","information_gain","importance","variance","mrmr","cmim","boruta"),
 genes = 50L,
 rf_selection_rule = "mean", # matches Appendix G; alternative: worst_then_mean
 # New orchestration preserves supplied fixed learner settings.
 # The historical nCV tuning archives are required for exact tuned-model reproduction.
 benchmark_mode = "supplied_fixed",
 outer_folds = 5L,
 candidates = expand.grid(
   representation=c("VST_fSVA","VST_ComBatRef","ComBatSeq_log1p"),
   filter=c("boruta","mrmr","cmim","importance","information_gain"),
   learner=c("Naive_Bayes","Linear_SVM","RBF_SVM","C5_0","Random_Forest","LightGBM","XGBoost"),
   genes=50L,stringsAsFactors=FALSE),
 # Optional list of candidate parameter lists per learner for reconstructed nested CV.
 # Supply actual archived candidates; no historical ranges are invented here.
 tuning_candidates = NULL,
 inner_folds = 3L,
 bootstrap_replicates = 2000L,
 calibration_bins = 6L, # RF: 6; matched architecture comparison: 8
 prediction_file = "results/TB_run_001/selected_holdout_predictions.csv"
)
