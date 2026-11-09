# Run from repository root: Rscript tests/test_core.R
source("R/core.R")
trace_rng();a<-runif(5);trace_rng();stopifnot(identical(a,runif(5)))
stopifnot(trace_seed(1,2)==as.integer((7832+1009+2*2018)%%2000000000))
set.seed(1)
x<-matrix(as.integer(rpois(120*40,20)),120,40,dimnames=list(paste0("g",1:120),paste0("s",1:40)))
m<-data.frame(sample_id=colnames(x),participant_id=paste0("p",1:40),study=rep(c("A","B"),each=20),disease_group="TB",disease=rep(c("case","control"),20))
m$outcome <- factor(m$disease,levels = c("control", "case"))
f<-trace_folds(m,5);stopifnot(length(f)==40,all(tabulate(f,5)==8))
z<-trace_representation("logCPM",x[,1:30],x[,31:40],m[1:30,],m[31:40,])
stopifnot(max(abs(z$new-log2(sweep(x[,31:40],2,colSums(x[,31:40]),"/")*1e6+1)))<1e-12)
y<-x[,31:40];y[,1]<-y[,1]*20
z2<-trace_representation("logCPM",x[,1:30],y,m[1:30,],m[31:40,])
stopifnot(identical(z$train,z2$train),identical(z$new[,-1],z2$new[,-1]))
if(requireNamespace("pROC",quietly=TRUE)&&requireNamespace("PRROC",quietly=TRUE)) {
 q<-trace_metrics(c("control","case"),c(.1,.9));stopifnot(q["AUROC"]==1,q["bACC"]==1,abs(q["Brier"]-.01)<1e-12)
}
message("Core tests passed")
