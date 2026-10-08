source("R/core.R");source("R/study_models.R")
a<-commandArgs(trailingOnly=TRUE);d<-trace_load(if(length(a))a[1] else "config/local.R")
out<-trace_record(d,"05_study_transportability")
i<-which(d$metadata$partition=="development");h<-which(d$metadata$partition=="holdout")
m<-d$metadata[i,,drop=FALSE];x<-d$counts[,i,drop=FALSE]
grid<-expand.grid(representation=c("VST_only","logCPM"),genes=c(25L,50L,100L),weighted=c(FALSE,TRUE),stringsAsFactors=FALSE)
pp<-list();ss<-list()
for(c in seq_len(nrow(grid))) {
 for(s in unique(m$study)) {
   tr<-which(m$study!=s);va<-which(m$study==s)
   if(length(unique(m$disease[tr]))!=2L||length(unique(m$disease[va]))!=2L)stop("Both classes required in each source and target")
   z<-trace_rf(x[,tr,drop=FALSE],x[,va,drop=FALSE],m[tr,,drop=FALSE],m[va,,drop=FALSE],grid$representation[c],grid$genes[c],grid$weighted[c],trace_seed(c,match(s,unique(m$study)),base=d$config$seed))
   pp[[length(pp)+1]]<-data.frame(candidate=c,study=s,sample_id=m$sample_id[va],disease=m$disease[va],probability=z$probability)
   ss[[length(ss)+1]]<-data.frame(candidate=c,study=s,as.list(trace_metrics(m$disease[va],z$probability)))
   saveRDS(z$genes,file.path(out,paste0("candidate_",c,"_study_",match(s,unique(m$study)),"_genes.rds")))
 }
}
p<-do.call(rbind,pp);s<-do.call(rbind,ss)
rank<-do.call(rbind,lapply(split(s,s$candidate),function(z)data.frame(candidate=z$candidate[1],mean_AUROC=mean(z$AUROC),worst_AUROC=min(z$AUROC),SD=sd(z$AUROC),mean_bACC=mean(z$bACC))))
# Explicit choice resolves the documented discrepancy without rewriting history.
if(is.null(d$config$rf_selection_rule))stop("Set rf_selection_rule to mean or worst_then_mean after reviewing candidate results")
write.csv(s,file.path(out,"study_metrics.csv"),row.names=FALSE)
write.csv(p,file.path(out,"study_predictions.csv"),row.names=FALSE)
if(d$config$rf_selection_rule=="mean")rank<-rank[order(-rank$mean_AUROC,rank$SD,rank$candidate),] else if(d$config$rf_selection_rule=="worst_then_mean")rank<-rank[order(-rank$worst_AUROC,-rank$mean_AUROC,rank$SD,-rank$mean_bACC),] else stop("Unknown RF selection rule")
best<-rank$candidate[1];g<-grid[best,];threshold<-trace_threshold(p[p$candidate==best,])
z<-trace_rf(x,d$counts[,h,drop=FALSE],m,d$metadata[h,,drop=FALSE],g$representation,g$genes,g$weighted,trace_seed(500,best,base=d$config$seed))
write.csv(rank,file.path(out,"candidate_ranking.csv"),row.names=FALSE)
write.csv(data.frame(sample_id=d$metadata$sample_id[h],disease=d$metadata$disease[h],probability=z$probability,threshold=threshold),file.path(out,"holdout_predictions.csv"),row.names=FALSE)
saveRDS(list(spec=g,threshold=threshold,fit=z),file.path(out,"final_rf.rds"))
