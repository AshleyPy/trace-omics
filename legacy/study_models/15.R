# ======================================================================
# PART 13
# DALEX EXPLAINER
# ======================================================================

dalex_x <- as.data.frame(
  final_x_ice,
  check.names = FALSE
)


# Ensure ICE rows are in metadata order if sample row names are available.
if (
  !is.null(
    rownames(
      dalex_x
    )
  ) &&
    all(
      final_ice_metadata$sample_id %in%
        rownames(
          dalex_x
        )
    )
) {

  dalex_x <- dalex_x[
    final_ice_metadata$sample_id,
    ,
    drop = FALSE
  ]
}


dalex_y <- as.integer(
  final_ice_metadata$outcome ==
    POSITIVE_CLASS
)


stopifnot(
  nrow(
    dalex_x
  ) ==
    length(
      dalex_y
    )
)


stopifnot(
  identical(
    colnames(
      dalex_x
    ),
    FINAL_GENES
  )
)


dalex_predict_function <- function(
  model,
  newdata
) {

  probability_matrix <- predict(
    model,
    data = as.data.frame(
      newdata,
      check.names = FALSE
    )
  )$predictions

  as.numeric(
    probability_matrix[
      ,
      POSITIVE_CLASS
    ]
  )
}


tb_rf_explainer <- DALEX::explain(
  model = FINAL_RF,
  data = dalex_x,
  y = dalex_y,
  predict_function = dalex_predict_function,
  label = paste0(
    "TB logCPM RF | ",
    LOCKED_TOP_N,
    " genes"
  ),
  type = "classification",
  verbose = FALSE
)


saveRDS(
  tb_rf_explainer,
  file.path(
    DALEX_DIR,
    "TB_final_RF_DALEX_explainer.rds"
  )
)
