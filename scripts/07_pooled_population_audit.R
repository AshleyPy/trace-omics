# Audit explicit Model A/B/shared-holdout manifests; never infer overlap from accessions.
source("R/core.R")
a<-commandArgs(trailingOnly=TRUE);d<-trace_load(if(length(a))a[1] else "config/local.R")
out<-trace_record(d,"07_pooled_population_audit");m<-d$metadata
if(!all(c("TB","RA","UC","CD","T1D")%in%m$disease_group))stop("All five disease labels required")
if(nrow(d$counts)!=4450L)stop("Main pooled comparison requires the frozen 4,450-gene universe")
write.csv(m[m$partition=="development" & m$disease_group!="TB",],file.path(out,"model_A_development.csv"),row.names=FALSE)
write.csv(m[m$partition=="development",],file.path(out,"model_B_development.csv"),row.names=FALSE)
write.csv(m[m$partition=="holdout",],file.path(out,"shared_holdout.csv"),row.names=FALSE)
for(group in unique(m$disease_group)) {
 source<-m$partition=="development" & m$disease_group!=group
 target<-m$partition=="development" & m$disease_group==group
 if(any(m$participant_id[source]%in%m$participant_id[target]))stop("Participant overlap across disease exclusion")
 write.csv(data.frame(sample_id=m$sample_id[source|target],role=ifelse(target[source|target],"target","source")),
   file.path(out,paste0("LODO_",group,".csv")),row.names=FALSE)
}
# Main VST+fSVA+mRMR production model is not reconstructed as an historical run.
# Use restored preprocessing states, selected genes and classifier objects to replay.
# The supplied alternative ComBat-Ref/Boruta/LightGBM implementation is in legacy/pooled_combatref.
