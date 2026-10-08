# ======================================================================
# PART 1
# ONE-TIME PYTHON ENVIRONMENT SETUP
#
# Run this ONCE.
#
# Requires Python >= 3.10.
# ======================================================================

if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

library(reticulate)


TB_PYTHON_ENV <- "tb_foundation_models"


# ======================================================================
# CREATE ENVIRONMENT IF IT DOES NOT EXIST
# ======================================================================

if (
  !TB_PYTHON_ENV %in%
    reticulate::virtualenv_list()
) {

  reticulate::virtualenv_create(
    envname = TB_PYTHON_ENV
  )
}


# ======================================================================
# ACTIVATE
# ======================================================================

reticulate::use_virtualenv(
  TB_PYTHON_ENV,
  required = TRUE
)


# ======================================================================
# INSTALL PACKAGES
# ======================================================================

reticulate::py_install(
  packages = c(
    "numpy",
    "scikit-learn",
    "lightgbm",
    "tabpfn"
  ),
  envname = TB_PYTHON_ENV,
  pip = TRUE
)


# ======================================================================
# AUDIT
# ======================================================================

print(
  reticulate::py_config()
)
