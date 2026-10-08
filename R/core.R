# New orchestration for TRACE-Omics; see docs/REPRODUCTION_STATUS.md.
# Functions operate on genes-by-samples count matrices and explicit manifests.
trace_seed <- function(..., base = 7832L) {
  z <- as.numeric(unlist(list(...)))
  as.integer((base + sum(z * (1009 * seq_along(z)))) %% 2000000000)
}
trace_rng <- function(seed = 7832L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
}
trace_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required packages: ", paste(absent, collapse = ", "))
}
trace_load <- function(config_file) {
  if (!file.exists(config_file)) stop("Configuration missing: ", config_file)
  cfg <- dget(config_file)
  required <- c("counts", "metadata", "partition", "output", "seed")
  if (!all(required %in% names(cfg))) stop("Configuration requires: ", paste(required, collapse=", "))
  for (key in c("counts", "metadata", "partition")) {
    if (!file.exists(cfg[[key]])) stop("Missing ", key, ": ", cfg[[key]])
  }
  x <- readRDS(cfg$counts)
  m <- read.csv(cfg$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "participant_id", "study", "disease_group", "disease")
  if (!all(required %in% names(m))) stop("Metadata requires: ", paste(required, collapse=", "))
  if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) || is.null(colnames(x)))
    stop("Counts must be a named numeric matrix, genes by samples")
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(abs(x - round(x)) > 1e-8))
    stop("Counts must contain finite non-negative integers")
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) || anyDuplicated(m$sample_id))
    stop("Duplicate gene or sample identifiers require resolution before modelling")
  if (anyNA(m[,required]) || any(vapply(m[,required], function(z) any(!nzchar(trimws(z))), logical(1))))
    stop("Metadata identifiers and labels must be populated")
  if (!setequal(colnames(x),m$sample_id)) stop("Matrix and metadata sample sets differ")
  m <- m[match(colnames(x),m$sample_id),,drop=FALSE]
  if (!all(m$disease %in% c("control","case"))) stop("Use exact disease labels control and case")
  m$outcome <- factor(m$disease, levels=c("control","case"))
  m$cohort_group <- factor(m$study)
  m$cohort <- factor(m$study) # adapter for supplied functions; scientific term is study
  if (any(colSums(x) == 0)) stop("Zero library size")
  z <- read.csv(cfg$partition, stringsAsFactors=FALSE)
  if (!all(c("sample_id","partition") %in% names(z)) || anyDuplicated(z$sample_id) ||
      !setequal(z$sample_id,m$sample_id)) stop("Partition manifest must identify every sample exactly once")
  m$partition <- z$partition[match(m$sample_id,z$sample_id)]
  if (!all(m$partition %in% c("development","holdout"))) stop("Partition labels: development, holdout")
  for (id in unique(m$participant_id)) {
    if (length(unique(m$partition[m$participant_id == id])) != 1L)
      stop("Participant crosses development/holdout: ",id)
  }
  if (anyDuplicated(m$participant_id))
    stop("Resolve repeated participant records explicitly before this runner; see data contract")
  for (part in c("development","holdout")) {
    if (length(unique(m$disease[m$partition == part])) != 2L) stop("Both classes required in ",part)
  }
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  trace_rng(cfg$seed)
  list(config=cfg,counts=x,metadata=m)
}
trace_record <- function(d,stage) {
  target <- file.path(d$config$output,stage)
  if (dir.exists(target)) stop("Stage output already exists: ",target,". Use a new output path for a new run.")
  dir.create(target,recursive=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(target,"sessionInfo.txt"))
  dput(d$config,file.path(target,"config.R"))
  saveRDS(list(RNGkind=RNGkind(),seed=d$config$seed,state=.Random.seed),file.path(target,"rng.rds"))
  write.csv(data.frame(path=unlist(d$config[c("counts","metadata","partition")]),
                      md5=unname(tools::md5sum(unlist(d$config[c("counts","metadata","partition")])))),
            file.path(target,"input_checksums.csv"),row.names=FALSE)
  target
}
trace_source_functions <- function() {
  e <- new.env(parent=globalenv())
  e$SEED <- 7832L; e$POSITIVE_CLASS <- "case"; e$NEGATIVE_CLASS <- "control"
  e$MIN_TOTAL_TRAINING_COUNT <- 10L; e$TOP_N <- 50L; e$MAX_SV <- 50L
  sys.source("R/supplied_functions.R",envir=e)
  e
}
trace_representation <- function(id, train, new=NULL, train_meta, new_meta=NULL, seed=7832L) {
  e <- trace_source_functions(); e$SEED <- seed
  if (id == "RAW_counts") return(list(train=train,new=new,state=list(genes=rownames(train))))
  if (id == "logCPM") {
    keep <- rowSums(train) >= 10
    transform <- function(x) log2(sweep(x[keep,,drop=FALSE],2,colSums(x),"/") * 1e6 + 1)
    if (any(colSums(train[keep,,drop=FALSE])==0) || (!is.null(new) && any(colSums(new[keep,,drop=FALSE])==0)))
      stop("Zero library size after training gene filter")
    return(list(train=transform(train),new=if(!is.null(new)) transform(new),state=list(genes=rownames(train)[keep])))
  }
  trace_require(c("DESeq2","SummarizedExperiment","matrixStats"))
  if (id == "VST_only") return(e$fit_apply_frozen_vst(train,new,train_meta))
  if (id == "VST_ComBat") {
    if (!is.null(new)) stop("Standard ComBat is descriptive-only in this driver")
    trace_require("neuroCombat")
    v <- e$fit_apply_frozen_vst(train,NULL,train_meta)
    fit <- neuroCombat::neuroCombat(dat=v$train,batch=train_meta$cohort,
        mod=NULL,parametric=TRUE)
    return(list(train=fit$dat.combat,new=NULL,state=fit$estimates))
  }
  if (id == "VST_fSVA") trace_require("sva")
  if (id == "VST_ComBatRef") {
    trace_require("neuroCombat")
    if (!is.null(new_meta) && any(!new_meta$study %in% train_meta$study))
      stop("ComBat-Ref requires evaluation batches represented in training")
  }
  if (id == "ComBatSeq_log1p") {
    trace_require(c("edgeR","matrixStats"))
    if (!is.null(new_meta) && any(!new_meta$study %in% train_meta$study))
      stop("Frozen ComBat-Seq requires represented batches")
  }
  e$fit_apply_representation(id,train,new,train_meta,new_meta,seed)
}
trace_select <- function(train,new,meta,method,k,seed) {
  trace_require(c("mlr3","mlr3filters","mlr3learners","matrixStats","data.table"))
  e <- trace_source_functions(); e$TOP_N <- k
  e$select_top_features(train,new,meta,method,k,seed)
}
trace_folds <- function(meta,k=5L,seed=7832L) {
  # New fold construction. Use saved historical folds for numerical reproduction.
  trace_rng(seed)
  if (anyDuplicated(meta$participant_id)) stop("Repeated participants require grouped resampling")
  strata <- interaction(meta$study,meta$disease,drop=TRUE)
  f <- integer(nrow(meta))
  for (idx in split(seq_len(nrow(meta)),strata)) {
    ids <- idx[sample.int(length(idx))]
    assignment <- rep(sample.int(k),length.out=length(ids))
    f[ids] <- assignment
  }
  for (j in seq_len(k)) {
    if (length(unique(meta$disease[f==j]))!=2L || length(unique(meta$disease[f!=j]))!=2L)
      stop("A fold lacks a class; provide a valid saved fold manifest or reduce folds")
  }
  f
}
trace_auc <- function(y,p) {
  trace_require("pROC")
  if (length(unique(y))<2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y,p,levels=c("control","case"),direction="<",quiet=TRUE)))
}
trace_metrics <- function(y,p,threshold=.5) {
  trace_require(c("pROC","PRROC"))
  if (length(y)!=length(p) || anyNA(p) || any(!is.finite(p)) || any(p<0|p>1)) stop("Invalid predictions")
  truth <- as.character(y)=="case"; positive <- p >= threshold
  sens <- if(any(truth)) mean(positive[truth]) else NA_real_
  spec <- if(any(!truth)) mean(!positive[!truth]) else NA_real_
  tp <- sum(truth & positive); fp <- sum(!truth & positive); fn <- sum(truth & !positive)
  c(AUROC=trace_auc(y,p), AUPRC=if(length(unique(truth))==2L)
      PRROC::pr.curve(scores.class0=p[truth],scores.class1=p[!truth],curve=FALSE)$auc.integral else NA_real_,
    bACC=mean(c(sens,spec)),Sens=sens,Spec=spec,F1=if(2*tp+fp+fn>0)2*tp/(2*tp+fp+fn) else NA_real_,
    Brier=mean((as.numeric(truth)-p)^2), prevalence=mean(truth))
}
trace_intervals <- function(y,p,threshold=.5,B=2000L,seed=7832L) {
  point <- trace_metrics(y,p,threshold); point <- point[names(point)!="prevalence"]; trace_rng(seed)
  if(length(unique(y))!=2L) stop("Both classes required for confidence intervals")
  r <- pROC::roc(y,p,levels=c("control","case"),direction="<",quiet=TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(r,method="delong"))
  classes <- split(seq_along(y),y)
  boot <- replicate(B,{
    ii <- unlist(lapply(classes,function(z)z[sample.int(length(z),length(z),replace=TRUE)]),use.names=FALSE)
    trace_metrics(y[ii],p[ii],threshold)[names(point)]
  })
  q <- t(apply(boot,1,quantile,probs=c(.025,.975),na.rm=TRUE,names=FALSE))
  q["AUROC",] <- auc_ci[c(1,3)]
  data.frame(metric=names(point),estimate=unname(point),lower=q[,1],upper=q[,2],
    method=ifelse(names(point)=="AUROC","DeLong","class-stratified percentile bootstrap"),row.names=NULL)
}
trace_calibration <- function(y,p,bins=8L) {
  # Equal-width quantiles with tied probabilities kept together.
  breaks <- unique(quantile(p,seq(0,1,length.out=bins+1),names=FALSE))
  group <- if(length(breaks)<2) rep(1L,length(p)) else cut(p,breaks,include.lowest=TRUE,labels=FALSE)
  do.call(rbind,lapply(split(seq_along(p),group),function(i){
    cases <- sum(y[i]=="case"); ci <- binom.test(cases,length(i))$conf.int
    data.frame(n=length(i),predicted=mean(p[i]),observed=cases/length(i),lower=ci[1],upper=ci[2])
  }))
}
trace_predict_pipeline <- function(train,new,meta,new_meta,spec,seed) {
  e <- trace_source_functions()
  r <- trace_representation(spec$representation,train,new,meta,new_meta,seed)
  f <- trace_select(r$train,r$new,meta,spec$filter,as.integer(spec$genes),seed)
  learner <- e$make_learner(spec$learner)
  if (!is.null(spec$parameters)) learner$param_set$values <- modifyList(learner$param_set$values,spec$parameters)
  task <- e$make_classification_task(f$train,meta,"trace_training")
  set.seed(seed); learner$train(task)
  x <- as.data.frame(t(f$validation));names(x)<-make.names(rownames(f$validation),unique=TRUE)
  probability <- learner$predict_newdata(x)$prob[,"case"]
  list(probability=probability,genes=f$genes,learner=learner,representation_state=r$state)
}
