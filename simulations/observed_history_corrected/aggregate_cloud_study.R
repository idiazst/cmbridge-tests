library(data.table)
outdir<-Sys.getenv('SIM_OUTPUT');incoming<-Sys.getenv('CLOUD_INCOMING')
dir.create(file.path(outdir,'jobs'),recursive=TRUE,showWarnings=FALSE)
dirs<-list.dirs(incoming,recursive=FALSE,full.names=TRUE)
dirs<-dirs[file.exists(file.path(dirs,'design.csv'))]
if(!length(dirs))stop('No shard checkpoints were uploaded.')
fingerprints<-vapply(dirs,function(d)readLines(file.path(d,'source-fingerprint.txt'))[1L],'')
stopifnot(length(unique(fingerprints))==1L)
metadata<-c('design.csv','specifications.csv','design-details.rds','source-fingerprint.txt')
for(file in metadata)stopifnot(file.copy(file.path(dirs[1L],file),file.path(outdir,file),overwrite=TRUE))
design<-fread(file.path(outdir,'design.csv'))
stopifnot(nrow(design)==1200L,uniqueN(design$job)==1200L,all(design[,.N,by=.(mechanism,n)]$N==200L))
seen<-character()
for(d in dirs) {
  for(file in list.files(d,'^(session|started|finished)-shard-',full.names=TRUE))
    stopifnot(file.copy(file,outdir,overwrite=TRUE))
  for(path in list.files(file.path(d,'jobs'),'\\.rds$',full.names=TRUE)) {
    job<-sub('\\.rds$','',basename(path));stopifnot(!job%in%seen,job%in%design$job)
    seen<-c(seen,job);x<-readRDS(path);wanted_job<-job;config<-design[job==wanted_job]
    stopifnot(identical(x$fingerprint,fingerprints[1L]),x$design$seed==config$seed,
      x$design$n==config$n,x$design$mechanism==config$mechanism,x$design$replicate==config$replicate)
    stopifnot(file.copy(path,file.path(outdir,'jobs',basename(path)),overwrite=FALSE))
  }
}
writeLines(c(paste('Saved',length(seen),'of 1200 planned datasets'),
  paste('Missing shards:',paste(setdiff(1:40,as.integer(sub('.*-','',basename(dirs)))),collapse=', ')),
  'Failures remain in checkpoints and are included in the summaries.'),file.path(outdir,'PROGRESS.txt'))
writeLines(if(length(seen)==1200L)'All planned datasets have saved checkpoints.'else
  'Partial study. Incomplete replications are not final results.',file.path(outdir,'STATUS.txt'))
if(length(seen)==1200L)writeLines(format(Sys.time(),tz='UTC',usetz=TRUE),file.path(outdir,'run-finished-at.txt'))
for(script in c('summarize.R','write_reduced_report.R'))
  source(file.path(Sys.getenv('STUDY_SOURCE'),script),local=new.env(parent=globalenv()))
cat(length(seen),'of 1200 checkpoints summarized.\n')
