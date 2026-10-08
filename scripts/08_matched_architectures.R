# New orchestration of the supplied architecture functions and ten candidate settings.
# Requires an already configured, version-matched Python environment and checkpoint.
source("R/core.R")
a<-commandArgs(trailingOnly=TRUE);d<-trace_load(if(length(a))a[1] else "config/local.R")
trace_require(c("reticulate","data.table","ranger","matrixStats","pROC","PRROC"))
if(is.null(d$config$python)||!file.exists(d$config$python))stop("Set python to the existing Python 3.10 executable")
reticulate::use_python(d$config$python,required=TRUE)
e<-new.env(parent=globalenv());e$POSITIVE_CLASS<-"case";e$NEGATIVE_CLASS<-"control"
e$RF_PREFILTER<-3000L;e$RF_IMPORTANCE_TREES<-250L;e$TABPFN_DEVICE<-"cpu"
e$np<-reticulate::import("numpy",convert=FALSE)
e$lightgbm_py<-reticulate::import("lightgbm",convert=FALSE)
e$tabpfn_py<-reticulate::import("tabpfn",convert=FALSE)
versions<-reticulate::py_run_string("import sys\nimport importlib.metadata as md\nversions = {p: md.version(p) for p in ['numpy','scikit-learn','lightgbm','torch','tabpfn']}\npython_version = sys.version",convert=TRUE)
expected<-c(numpy="1.22.4",`scikit-learn`="1.2.2",lightgbm="4.6.0",torch="2.4.1",tabpfn="0.1.11")
actual<-unlist(versions$versions)
if(any(sub("[+].*$","",actual[names(expected)])!=expected))stop("Python package versions differ from the recorded specification")
if(!startsWith(versions$python_version,"3.10."))stop("Use recorded Python 3.10")
sys.source("R/supplied_architecture_functions.R",envir=e)
out<-trace_record(d,"08_matched_architectures")
writeLines(capture.output(versions$versions),file.path(out,"python_versions.txt"))
i<-which(d$metadata$partition=="development");h<-which(d$metadata$partition=="holdout")
m<-d$metadata[i,,drop=FALSE];x<-d$counts[,i,drop=FALSE]
prepare<-function(tr,va,xx=x,mm=m){
 z<-trace_representation("logCPM",xx[,tr,drop=FALSE],xx[,va,drop=FALSE],mm[tr,,drop=FALSE],mm[va,,drop=FALSE])
 genes<-e$select_foundation_rf_features(z$train,mm[tr,,drop=FALSE],100L,trace_seed(length(tr),length(va),base=d$config$seed))
 list(train=t(z$train[genes,,drop=FALSE]),new=t(z$new[genes,,drop=FALSE]),genes=genes)
}
tune<-function(tr,seed){
 f<-trace_folds(m[tr,,drop=FALSE],3L,seed);score<-matrix(NA_real_,nrow(e$lightgbm_grid),3L)
 for(g in 1:3){ti<-tr[f!=g];vi<-tr[f==g];z<-prepare(ti,vi)
  for(k in seq_len(nrow(e$lightgbm_grid))){fit<-e$fit_lightgbm_object(z$train,m$outcome[ti],e$lightgbm_grid[k],trace_seed(seed,g,k));p<-e$predict_python_classifier(fit,z$new);score[k,g]<-trace_auc(m$disease[vi],p)}
 }
 list(index=which.max(rowMeans(score)),scores=score)
}
run_pair<-function(tr,va,label,xx=x,mm=m,selection=NULL){
 if(is.null(selection))selection<-tune(tr,trace_seed(length(tr),base=d$config$seed))
 z<-prepare(tr,va,xx,mm);seed<-trace_seed(length(tr),length(va),base=d$config$seed)
 fits<-list(LightGBM=e$fit_lightgbm_object(z$train,mm$outcome[tr],e$lightgbm_grid[selection$index],seed),TabPFN=e$fit_tabpfn_object(z$train,mm$outcome[tr],seed))
 saveRDS(list(genes=z$genes,selection=selection),file.path(out,paste0(label,"_specification.rds")))
 do.call(rbind,lapply(names(fits),function(name){
  p<-e$predict_python_classifier(fits[[name]],z$new)
  if(label=="holdout")reticulate::py_save_object(fits[[name]],file.path(out,paste0(name,".pkl")))
  data.frame(evaluation=label,model=name,sample_id=mm$sample_id[va],disease=mm$disease[va],probability=p,threshold=.5)
 }))
}
fold<-trace_folds(m,5L,d$config$seed);predictions<-list()
write.csv(data.frame(sample_id=m$sample_id,fold=fold),file.path(out,"outer_folds.csv"),row.names=FALSE)
for(f in 1:5)predictions[[length(predictions)+1]]<-run_pair(which(fold!=f),which(fold==f),paste0("outer_",f))
for(s in unique(m$study))predictions[[length(predictions)+1]]<-run_pair(which(m$study!=s),which(m$study==s),paste0("study_",match(s,unique(m$study))))
selection<-tune(seq_len(nrow(m)),trace_seed(999,base=d$config$seed))
# Shared full-data ordering retains the original development/holdout manifest.
predictions[[length(predictions)+1]]<-run_pair(i,h,"holdout",d$counts,d$metadata,selection)
write.csv(do.call(rbind,predictions),file.path(out,"all_predictions.csv"),row.names=FALSE)
