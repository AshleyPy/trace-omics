# ======================================================================
# PART 15
# REPEATED Delta-AUROC PERMUTATION IMPORTANCE
#
# Delta AUROC =
#
#     baseline AUROC
#       -
#     AUROC after permuting one gene
#
# Positive:
#   gene supports discrimination.
#
# Zero:
#   little unique contribution.
#
# Negative:
#   permutation improved discrimination.
#
# This is POST-HOC interpretation only.
# ======================================================================


baseline_probability <- dalex_predict_function(
  FINAL_RF,
  dalex_x
)


baseline_truth <- factor(
  ifelse(
    dalex_y == 1,
    POSITIVE_CLASS,
    NEGATIVE_CLASS
  ),
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  )
)


baseline_roc <- pROC::roc(
  response = baseline_truth,
  predictor = baseline_probability,
  levels = c(
    NEGATIVE_CLASS,
    POSITIVE_CLASS
  ),
  direction = "<",
  quiet = TRUE
)


baseline_auc <- as.numeric(
  pROC::auc(
    baseline_roc
  )
)


cat(
  "\nBaseline diagnostic ICE AUROC: ",
  sprintf(
    "%.4f",
    baseline_auc
  ),
  "\n",
  sep = ""
)


importance_raw_rows <- vector(
  mode = "list",
  length = ncol(
    dalex_x
  )
)


for (
  gene_index in seq_len(
    ncol(
      dalex_x
    )
  )
) {

  current_gene <- colnames(
    dalex_x
  )[[gene_index]]

  original_values <- dalex_x[
    ,
    gene_index
  ]

  gene_rows <- vector(
    mode = "list",
    length = B_IMPORTANCE
  )

  for (
    permutation_index in seq_len(
      B_IMPORTANCE
    )
  ) {

    set.seed(
      SEED_PUBLICATION +
        gene_index *
          10000L +
        permutation_index
    )

    x_permuted <- dalex_x

    x_permuted[
      ,
      gene_index
    ] <- sample(
      original_values,
      size = length(
        original_values
      ),
      replace = FALSE
    )

    permuted_probability <- dalex_predict_function(
      FINAL_RF,
      x_permuted
    )

    permuted_roc <- pROC::roc(
      response = baseline_truth,
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

    gene_rows[[permutation_index]] <-
      data.table::data.table(

        Gene =
          current_gene,

        Permutation =
          permutation_index,

        Baseline_AUROC =
          baseline_auc,

        Permuted_AUROC =
          permuted_auc,

        Delta_AUROC =
          baseline_auc -
          permuted_auc
      )
  }

  importance_raw_rows[[gene_index]] <-
    data.table::rbindlist(
      gene_rows
    )

  if (
    gene_index %% 5L == 0L ||
      gene_index ==
        ncol(
          dalex_x
        )
  ) {

    message(
      "Permutation importance: ",
      gene_index,
      "/",
      ncol(
        dalex_x
      ),
      " genes complete."
    )
  }
}


importance_raw_dt <- data.table::rbindlist(
  importance_raw_rows,
  use.names = TRUE,
  fill = TRUE
)


importance_summary_dt <- importance_raw_dt[
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
        Delta_AUROC > 0,
        na.rm = TRUE
      )
  ),
  by = Gene
]


data.table::setorder(
  importance_summary_dt,
  -Mean_Delta_AUROC,
  Gene
)


importance_summary_dt[
  ,
  Rank := seq_len(
    .N
  )
]


print(
  importance_summary_dt[
    1:min(
      20L,
      .N
    )
  ]
)


data.table::fwrite(
  importance_raw_dt,
  file.path(
    DALEX_DIR,
    "permutation_importance_raw.csv"
  )
)


data.table::fwrite(
  importance_summary_dt,
  file.path(
    TABLE_DIR,
    "Table_5_predictive_gene_importance.csv"
  )
)
