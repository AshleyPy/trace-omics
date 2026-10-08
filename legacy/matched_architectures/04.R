# ======================================================================
# RA TABPFN PERMUTATION IMPORTANCE
# SAFE CHECKPOINTED VERSION
#
# RUN IN THE CURRENT R SESSION
#
# DOES NOT RERUN:
#   - LightGBM importance
#   - preprocessing
#   - feature selection
#   - fitting
#   - diagnostic ICE
#
# Saves after EVERY completed gene.
# ======================================================================


# ======================================================================
# 1. SETTINGS
# ======================================================================

TABPFN_PERMUTATIONS <- 100L

SEED_TABPFN_IMPORTANCE <- SEED_FOUNDATION + 88200002L


TABPFN_CHECKPOINT_FILE <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_checkpoint.rds"
)


TABPFN_RAW_FILE <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_raw.csv"
)


TABPFN_SUMMARY_FILE <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_summary.csv"
)


TABPFN_FINAL_RDS <- file.path(
  DALEX_COMPARE_DIR,
  "RA_TabPFN_permutation_importance_complete.rds"
)


# ======================================================================
# 2. SUPPRESS PYTORCH WARNING NOISE
# ======================================================================

reticulate::py_run_string("
import warnings
warnings.filterwarnings('ignore')
")


# ======================================================================
# 3. ORIGINAL TABPFN BASELINE
# ======================================================================

truth_factor <- factor(
  ice_truth,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  )
)


tabpfn_baseline_probability <- predict_tabpfn_dalex(
  FINAL_TABPFN_MODEL,
  dalex_compare_x
)


tabpfn_baseline_roc <- pROC::roc(
  response = truth_factor,
  predictor = tabpfn_baseline_probability,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


TABPFN_BASELINE_AUROC <- as.numeric(
  pROC::auc(
    tabpfn_baseline_roc
  )
)


cat(
  "\n============================================================\n",
  "RA TABPFN CHECKPOINTED PERMUTATION IMPORTANCE\n",
  "============================================================\n",
  "Diagnostic ICE baseline AUROC: ",
  sprintf(
    "%.4f",
    TABPFN_BASELINE_AUROC
  ),
  "\nGenes: ",
  ncol(
    dalex_compare_x
  ),
  "\nPermutations per gene: ",
  TABPFN_PERMUTATIONS,
  "\n============================================================\n",
  sep = ""
)


# ======================================================================
# 4. LOAD EXISTING CHECKPOINT IF ONE EXISTS
# ======================================================================

if (
  file.exists(
    TABPFN_CHECKPOINT_FILE
  )
) {

  checkpoint_object <- readRDS(
    TABPFN_CHECKPOINT_FILE
  )


  tabpfn_raw_rows <- checkpoint_object$raw_rows


  completed_tabpfn_genes <- checkpoint_object$completed_genes


  message(
    "RESUME: ",
    length(
      completed_tabpfn_genes
    ),
    " genes already completed."
  )

} else {

  tabpfn_raw_rows <- list()


  completed_tabpfn_genes <- character()


  message(
    "Starting TabPFN importance from gene 1."
  )
}


# ======================================================================
# 5. GENES STILL TO RUN
# ======================================================================

all_tabpfn_genes <- colnames(
  dalex_compare_x
)


remaining_tabpfn_genes <- setdiff(
  all_tabpfn_genes,
  completed_tabpfn_genes
)


cat(
  "\nAlready complete: ",
  length(
    completed_tabpfn_genes
  ),
  "\nRemaining: ",
  length(
    remaining_tabpfn_genes
  ),
  "\n",
  sep = ""
)


# ======================================================================
# 6. RUN TABPFN IMPORTANCE
#
# IMPORTANT:
# Each permutation uses its own predict_proba() call.
# This preserves the original permutation-analysis behaviour exactly.
# ======================================================================

for (
  current_gene in remaining_tabpfn_genes
) {

  gene_index <- match(
    current_gene,
    all_tabpfn_genes
  )


  original_values <- dalex_compare_x[[current_gene]]


  current_gene_rows <- vector(
    mode = "list",
    length = TABPFN_PERMUTATIONS
  )


  cat(
    "\n------------------------------------------------------------\n",
    "TabPFN gene ",
    gene_index,
    "/",
    length(
      all_tabpfn_genes
    ),
    ": ",
    current_gene,
    "\n",
    "------------------------------------------------------------\n",
    sep = ""
  )


  for (
    permutation_index in seq_len(
      TABPFN_PERMUTATIONS
    )
  ) {

    set.seed(
      SEED_TABPFN_IMPORTANCE +
        gene_index *
          10000L +
        permutation_index
    )


    x_permuted <- dalex_compare_x


    x_permuted[[current_gene]] <- sample(
      original_values,
      size = length(
        original_values
      ),
      replace = FALSE
    )


    permuted_probability <- predict_tabpfn_dalex(
      FINAL_TABPFN_MODEL,
      x_permuted
    )


    permuted_roc <- pROC::roc(
      response = truth_factor,
      predictor = permuted_probability,
      levels = c(
        NEGATIVE_CLASS,
        POSITIVE_CLASS
      ),
      direction = "<",
      quiet = TRUE
    )


    permuted_auc <- as.numeric(
      pROC::auc(
        permuted_roc
      )
    )


    current_gene_rows[[permutation_index]] <-
      data.table::data.table(

        Model =
          "TabPFN",

        Gene =
          current_gene,

        Permutation =
          permutation_index,

        Baseline_AUROC =
          TABPFN_BASELINE_AUROC,

        Permuted_AUROC =
          permuted_auc,

        Delta_AUROC =
          TABPFN_BASELINE_AUROC -
          permuted_auc
      )


    if (
      permutation_index %% 20L ==
        0L
    ) {

      message(
        "  ",
        current_gene,
        ": ",
        permutation_index,
        "/",
        TABPFN_PERMUTATIONS,
        " permutations"
      )
    }
  }


  # ====================================================================
  # SAVE CURRENT GENE
  # ====================================================================

  tabpfn_raw_rows[[current_gene]] <-
    data.table::rbindlist(
      current_gene_rows
    )


  completed_tabpfn_genes <- c(
    completed_tabpfn_genes,
    current_gene
  )


  completed_tabpfn_genes <- unique(
    completed_tabpfn_genes
  )


  # --------------------------------------------------------------------
  # CHECKPOINT AFTER EVERY GENE
  # --------------------------------------------------------------------

  saveRDS(
    list(
      raw_rows =
        tabpfn_raw_rows,
      completed_genes =
        completed_tabpfn_genes,
      baseline_auroc =
        TABPFN_BASELINE_AUROC,
      permutations =
        TABPFN_PERMUTATIONS
    ),
    TABPFN_CHECKPOINT_FILE
  )


  # --------------------------------------------------------------------
  # ALSO WRITE CURRENT RAW CSV
  # --------------------------------------------------------------------

  current_raw_dt <- data.table::rbindlist(
    tabpfn_raw_rows,
    use.names = TRUE,
    fill = TRUE
  )


  data.table::fwrite(
    current_raw_dt,
    TABPFN_RAW_FILE
  )


  message(
    "SAVED: TabPFN ",
    length(
      completed_tabpfn_genes
    ),
    "/",
    length(
      all_tabpfn_genes
    ),
    " genes complete."
  )
}


# ======================================================================
# 7. COMBINE COMPLETED RAW RESULTS
# ======================================================================

tabpfn_importance_raw <- data.table::rbindlist(
  tabpfn_raw_rows,
  use.names = TRUE,
  fill = TRUE
)


# ======================================================================
# 8. SUMMARISE TABPFN IMPORTANCE
# ======================================================================

tabpfn_importance_summary <- tabpfn_importance_raw[
  ,
  .(

    Baseline_AUROC =
      unique(
        Baseline_AUROC
      )[[1L]],

    Mean_permuted_AUROC =
      mean(
        Permuted_AUROC,
        na.rm = TRUE
      ),

    Mean_Delta_AUROC =
      mean(
        Delta_AUROC,
        na.rm = TRUE
      ),

    SD_Delta_AUROC =
      stats::sd(
        Delta_AUROC,
        na.rm = TRUE
      ),

    Lower95_Delta_AUROC =
      as.numeric(
        stats::quantile(
          Delta_AUROC,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

    Upper95_Delta_AUROC =
      as.numeric(
        stats::quantile(
          Delta_AUROC,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

    Positive_permutation_fraction =
      mean(
        Delta_AUROC >
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    Model,
    Gene
  )
]


data.table::setorder(
  tabpfn_importance_summary,
  -Mean_Delta_AUROC,
  Gene
)


tabpfn_importance_summary[
  ,
  Rank :=
    seq_len(
      .N
    )
]


# ======================================================================
# 9. CREATE FINAL OBJECT
# ======================================================================

tabpfn_importance_result <- list(

  raw =
    tabpfn_importance_raw,

  summary =
    tabpfn_importance_summary
)


# ======================================================================
# 10. SAVE FINAL TABPFN RESULTS
# ======================================================================

saveRDS(
  tabpfn_importance_result,
  TABPFN_FINAL_RDS
)


data.table::fwrite(
  tabpfn_importance_raw,
  TABPFN_RAW_FILE
)


data.table::fwrite(
  tabpfn_importance_summary,
  TABPFN_SUMMARY_FILE
)


# ======================================================================
# 11. FINAL TABPFN AUDIT
# ======================================================================

stopifnot(
  nrow(
    tabpfn_importance_summary
  ) ==
    ncol(
      dalex_compare_x
    )
)


stopifnot(
  all(
    tabpfn_importance_raw[
      ,
      .N,
      by = Gene
    ]$N ==
      TABPFN_PERMUTATIONS
  )
)


cat(
  "\n\n============================================================\n",
  "PASS: RA TABPFN IMPORTANCE COMPLETE AND SAVED\n",
  "============================================================\n",
  "Genes: ",
  nrow(
    tabpfn_importance_summary
  ),
  "\nPermutations per gene: ",
  TABPFN_PERMUTATIONS,
  "\nTotal permutation predictions: ",
  nrow(
    tabpfn_importance_raw
  ),
  "\nBaseline diagnostic ICE AUROC: ",
  sprintf(
    "%.4f",
    TABPFN_BASELINE_AUROC
  ),
  "\n\nSaved summary:\n",
  TABPFN_SUMMARY_FILE,
  "\n============================================================\n",
  sep = ""
)
