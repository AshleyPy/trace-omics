files<-list.files(".",pattern="[.]R$",recursive=TRUE,full.names=TRUE)
for(f in files)tryCatch(parse(file=f),error=function(e)stop(f,": ",conditionMessage(e)))
message(length(files)," R files parsed")
