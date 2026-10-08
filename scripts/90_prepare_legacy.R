# Combine supplied chunks and apply explicitly reviewed, literal path replacements.
# Creates a script for inspection; does not execute it or change analytical settings.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=1L)stop("Usage: Rscript scripts/90_prepare_legacy.R config/legacy_local.R")
cfg<-dget(a[1]);dir<-file.path("legacy",cfg$workflow)
if(!dir.exists(dir))stop("Unknown workflow")
files<-sort(list.files(dir,pattern="[.]R$",full.names=TRUE))
if(!is.null(cfg$chunks))files<-files[cfg$chunks]
text<-paste(vapply(files,function(f)paste(readLines(f,warn=FALSE),collapse="\n"),character(1)),collapse="\n\n")
if(!is.null(cfg$replacements))for(k in names(cfg$replacements)) {
 if(!grepl(k,text,fixed=TRUE))stop("Replacement target missing: ",k)
 text<-gsub(k,cfg$replacements[[k]],text,fixed=TRUE)
}
if(file.exists(cfg$output))stop("Prepared output exists; choose a new filename")
dir.create(dirname(cfg$output),recursive=TRUE,showWarnings=FALSE)
writeLines(text,cfg$output)
parse(file=cfg$output)
dput(cfg,paste0(cfg$output,".preparation.R"))
message("Prepared for review: ",cfg$output,". Check data identity, path settings and package-installation sections before running.")
