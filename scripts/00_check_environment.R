# Audit only: no installation, environment changes or analysis execution.
cran <- c("data.table","matrixStats","cluster","mlr3","mlr3learners","mlr3extralearners",
 "mlr3filters","ranger","pROC","PRROC","Boruta","praznik","FSelectorRcpp",
 "e1071","naivebayes","C50","lightgbm","xgboost","reticulate","DALEX","ingredients","renv")
bioc <- c("DESeq2","SummarizedExperiment","sva","edgeR")
other <- "neuroCombat"
pkgs <- unique(c(cran,bioc,other))
versions <- vapply(pkgs,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else "MISSING",character(1))
print(data.frame(package=pkgs,version=unname(versions)),row.names=FALSE)
cat("\nR version:",R.version.string,"\n")
cat("Record sessionInfo() and create a real renv lock after verifying the working environment.\n")
cat("Use BiocManager for Bioconductor dependencies; see the neuroCombat R repository for its installation.\n")
cat("Install only the learner/filter dependencies used by your selected run.\n")
if(any(versions=="MISSING"))message("Some dependencies are missing; inspect the selected workflow before execution.")
