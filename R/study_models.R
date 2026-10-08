# New reusable RF orchestration with the recorded tree counts and feature counts.
trace_rf <- function(train,new,meta,new_meta,representation,genes,weighted,seed) {
 trace_require(c("ranger","matrixStats"))
 z <- trace_representation(representation,train,new,meta,new_meta,seed)
 v <- matrixStats::rowVars(z$train); keep <- which(is.finite(v)&v>0)
 keep <- keep[order(-v[keep],rownames(z$train)[keep])]
 keep <- head(keep,3000L)
 x <- as.data.frame(t(z$train[keep,,drop=FALSE])); colnames(x) <- make.names(colnames(x),unique=TRUE)
 if(ncol(x)<genes)stop("Fewer variable genes than requested")
 y <- factor(meta$disease,levels=c("control","case"))
 weights <- if(weighted) length(y)/(2*table(y)) else NULL
 screen <- ranger::ranger(x=x,y=y,probability=TRUE,num.trees=250L,
   mtry=max(1L,floor(sqrt(ncol(x)))),importance="permutation",class.weights=weights,num.threads=1L,seed=seed)
 imp <- screen$variable.importance
 ord <- order(-imp,names(imp),na.last=NA);take <- head(ord,genes)
 selected <- rownames(z$train)[keep[take]]
 xt <- x[,take,drop=FALSE]
 xv <- as.data.frame(t(z$new[selected,,drop=FALSE]));names(xv)<-names(xt)
 model <- ranger::ranger(x=xt,y=y,probability=TRUE,num.trees=1000L,
   mtry=max(1L,floor(sqrt(ncol(xt)))),min.node.size=5L,class.weights=weights,num.threads=1L,seed=seed)
 list(probability=as.numeric(predict(model,data=xv)$predictions[,"case"]),genes=selected,model=model,representation_state=z$state)
}
trace_threshold <- function(p) {
 grid <- seq(.05,.95,.01)
 scores <- t(vapply(grid,function(t){
   b <- vapply(split(p,p$study),function(z)trace_metrics(z$disease,z$probability,t)["bACC"],numeric(1))
   c(mean=mean(b),worst=min(b))
 },numeric(2)))
 grid[order(-scores[,1],-scores[,2],abs(grid-.5))[1]]
}
