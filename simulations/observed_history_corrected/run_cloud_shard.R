# Execution dispatcher only. Every fit is run by the unchanged frozen study.R.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
outdir <- Sys.getenv('SIM_OUTPUT')
frozen <- normalizePath(Sys.getenv('FROZEN_ROOT'))
dispatcher <- normalizePath(Sys.getenv('CLOUD_DISPATCHER'))
shards <- as.integer(Sys.getenv('SIM_SHARDS','40'))
shard <- as.integer(Sys.getenv('SIM_SHARD'))
workers <- as.integer(Sys.getenv('SIM_WORKERS','4'))
stopifnot(nzchar(outdir),shards==40L,shard>=1L,shard<=shards,workers==4L)
dir.create(file.path(outdir,'jobs'),recursive=TRUE,showWarnings=FALSE)
mechanisms <- c('binary_longitudinal','discrete_dose')
sizes <- c(500L,1000L,4000L)
grid <- CJ(mechanism=mechanisms,n=sizes,replicate=seq_len(200L))
grid[,seed:=5000000L+match(mechanism,mechanisms)*100000L+match(n,sizes)*1000L+replicate]
grid[,job:=sprintf('%s-n%d-r%03d',mechanism,n,replicate)]
ordered <- order(grid$replicate,grid$n,grid$mechanism)
grid[ordered,shard:=((seq_along(ordered)-1L)%%shards)+1L]
stopifnot(nrow(grid)==1200L,uniqueN(grid$job)==1200L,all(grid[,.N,by=.(mechanism,n)]$N==200L),
  all(grid[,.N,by=shard]$N==30L))
sources <- c(file.path(Sys.getenv('STUDY_SOURCE'),c('dgp.R','study.R','run_reduced.R')),
  list.files(file.path(frozen,'source/lmtp/R'),'\\.R$',full.names=TRUE),
  list.files(file.path(frozen,'source/cmbridge/R'),'\\.R$',full.names=TRUE),
  file.path(frozen,c('source/lmtp/DESCRIPTION','source/cmbridge/DESCRIPTION')))
hashes <- tools::md5sum(sort(sources));names(hashes)<-substring(names(hashes),nchar(frozen)+2L)
packages <- c('lmtp','cmbridge','data.table','SuperLearner','glmnet','earth','Matrix','digest')
versions <- vapply(packages,function(x)as.character(packageVersion(x)),'')
fingerprint <- digest::digest(list(sources=hashes,dispatcher=unname(tools::md5sum(dispatcher)),
  R=R.version$version.string,packages=versions),serializeVersion=2L)
manifest <- file.path(outdir,'source-fingerprint.txt')
if(file.exists(manifest)&&readLines(manifest)[1L]!=fingerprint)stop('Cannot mix checkpoints from different sources.')
writeLines(fingerprint,manifest)
fwrite(grid,file.path(outdir,'design.csv'))
fwrite(data.table(estimator=c('sdr','tmle'),beta_correct=TRUE,lambda_correct=TRUE,
  ratio_correct=TRUE,m_correct=TRUE),file.path(outdir,'specifications.csv'))
saveRDS(list(workers=80L,workers_per_shard=workers,shards=shards,replications=200L,sizes=sizes,
  mechanisms=mechanisms,dose_levels=0:3,all_correct=TRUE,learner_groups=3L,
  measurement_balanced=TRUE,visits='missing_both',execution='GitHub standard Ubuntu 24.04 runners'),
  file.path(outdir,'design-details.rds'))
writeLines(c(capture.output(sessionInfo()),paste('fingerprint',fingerprint),
  paste('shard',shard,'of',shards),paste('workers',workers)),file.path(outdir,paste0('session-shard-',shard,'.txt')))
started_path<-file.path(outdir,paste0('started-shard-',shard,'.txt'))
if(!file.exists(started_path))writeLines(format(Sys.time(),tz='UTC',usetz=TRUE),started_path)
selected <- ordered[grid$shard[ordered]==shard]
assigned <- length(selected)
limit <- as.integer(Sys.getenv('SIM_DATASET_LIMIT','30'))
stopifnot(limit>=1L,limit<=assigned)
selected <- head(selected,limit)
one <- function(i) {
  config<-grid[i];path<-file.path(outdir,'jobs',paste0(config$job,'.rds'))
  if(file.exists(path)) {
    previous<-readRDS(path)
    stopifnot(identical(previous$fingerprint,fingerprint),previous$design$seed==config$seed)
    return(TRUE)
  }
  cat('Starting',config$job,'worker',Sys.getpid(),'\n')
  started<-proc.time()[3L];warnings<-character()
  result<-tryCatch(withCallingHandlers(
    run_dataset(config$mechanism,config$n,config$replicate,config$seed,all_correct=TRUE,
      dose_max=3L,learner_groups=3L,balance_measurement=TRUE,visits='missing_both'),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e)list(design=as.list(config[,.(mechanism,n,replicate,seed)]),
      errors=data.table(error=conditionMessage(e))))
  result$fingerprint<-fingerprint;result$seconds<-proc.time()[3L]-started
  result$job_warnings<-unique(warnings)
  tmp<-paste0(path,'.tmp');saveRDS(result,tmp,compress='gzip');stopifnot(file.rename(tmp,path))
  status<-if(nrow(result$errors))paste(unique(result$errors$error),collapse=' | ')else 'OK'
  cat(config$job,sprintf('%.1fs',result$seconds),status,'\n');TRUE
}
for(batch in split(selected,ceiling(seq_along(selected)/(workers*3L)))) {
  answers<-parallel::mclapply(batch,one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  saved<-file.exists(file.path(outdir,'jobs',paste0(grid$job[batch],'.rds')))
  if(length(answers)!=length(batch)||!all(vapply(answers,isTRUE,FALSE))||!all(saved))
    stop('A worker did not save its checkpoint. Retain output and resume its unfinished jobs.')
  writeLines(c(paste('shard',shard),paste('completed',length(list.files(file.path(outdir,'jobs'),'\\.rds$')),
    'of',assigned),format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(outdir,'PROGRESS.txt'))
}
if(limit==assigned) {
  writeLines(format(Sys.time(),tz='UTC',usetz=TRUE),file.path(outdir,paste0('finished-shard-',shard,'.txt')))
  cat('All 30 assigned datasets saved. Final-estimator failures remain recorded in checkpoints.\n')
} else cat('Initial four-dataset checkpoint batch saved; the shard is incomplete.\n')
