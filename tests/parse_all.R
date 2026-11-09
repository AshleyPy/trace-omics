# Parse R source without loading packages or reading patient data.
# Executed from the repository root by GitHub Actions.
if (!file.exists("R/core.R")) stop("Run this check from TRACE-Omics project root.")
files <- sort(list.files(".", pattern = "[.]R$", recursive = TRUE, full.names = TRUE))
files <- files[!grepl("(^|/)(renv/library|renv/staging|[.]git|[.]Rproj.user)/", files)]
if (!length(files)) stop("No R scripts found.")
failures <- character()
for (file in files) {
  err <- tryCatch({
    parse(file = file, keep.source = FALSE)
    NULL
  }, error = function(e) conditionMessage(e))
  if (is.null(err)) {
    cat("PASS ", file, "\n", sep = "")
  } else {
    cat("FAIL ", file, ": ", err, "\n", sep = "")
    failures <- c(failures, file)
  }
}
cat("Parsed ", length(files), " R scripts; failures: ", length(failures), "\n", sep = "")
if (length(failures)) stop("Fix the listed R parser errors before release.", call. = FALSE)
