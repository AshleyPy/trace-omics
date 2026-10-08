# New allocation implementation; original sampling code is in legacy/equal_size/01.R.
# Historical comparisons require the ORIGINAL selected IDs, not a newly generated draw.
source("R/core.R")
a<-commandArgs(trailingOnly=TRUE);d<-trace_load(if(length(a))a[1] else "config/local.R")
out<-trace_record(d,"06_equal_size_manifest");m<-d$metadata
if(length(unique(m$disease_group))!=1L)stop("Run separately per disease group")
allocate<-function(sizes,total){
 if(sum(sizes)<total)stop("Insufficient class size")
 quota<-total*sizes/sum(sizes);n<-floor(quota)
 if(sum(n)<total){ii<-order(-(quota-n),names(sizes));n[ii[seq_len(total-sum(n))]]<-n[ii[seq_len(total-sum(n))]]+1L}
 n
}
selected<-integer(); holdout<-integer()
for(label in c("control","case")) {
 ii<-which(m$disease==label);strata<-split(ii,m$study[ii]);quota<-allocate(lengths(strata),41L)
 take<-unlist(Map(function(idx,n)idx[sample.int(length(idx),n)],strata,quota),use.names=FALSE)
 by_study<-split(take,m$study[take]);test_quota<-allocate(lengths(by_study),8L)
 test<-unlist(Map(function(idx,n)idx[sample.int(length(idx),n)],by_study,test_quota),use.names=FALSE)
 selected<-c(selected,take);holdout<-c(holdout,test)
}
mm<-m[selected,c("sample_id","participant_id","study","disease_group","disease")]
manifest<-data.frame(sample_id=mm$sample_id,partition=ifelse(selected%in%holdout,"holdout","development"))
write.csv(mm,file.path(out,"metadata.csv"),row.names=FALSE)
write.csv(manifest,file.path(out,"partition.csv"),row.names=FALSE)
saveRDS(d$counts[,mm$sample_id,drop=FALSE],file.path(out,"counts.rds"))
