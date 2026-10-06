# Instrument one frozen-data fit for runtime diagnosis. This is never pooled
# into the cloud study and does not change any fitting setting or sample split.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
outdir<-Sys.getenv('RUNTIME_PROFILE_OUTPUT');stopifnot(nzchar(outdir))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
mechanism<-Sys.getenv('PROFILE_MECHANISM','binary_longitudinal')
n<-as.integer(Sys.getenv('PROFILE_N','4000'))
replicate<-as.integer(Sys.getenv('PROFILE_REPLICATE','2'))
stopifnot(mechanism%in%c('binary_longitudinal','discrete_dose'),n%in%c(500L,1000L,4000L),
  replicate>=1L,replicate<=200L)
seed<-5000000L+match(mechanism,c('binary_longitudinal','discrete_dose'))*100000L+
  match(n,c(500L,1000L,4000L))*1000L+replicate
stopifnot(as.character(packageVersion('cmbridge'))=='0.3.0.9020',
  as.character(packageVersion('lmtp'))=='1.6.0.9021')
result_file<-file.path(outdir,'result.rds')
if(file.exists(result_file))stop('A saved runtime-profile result already exists; do not overwrite it.')
writeLines(c(capture.output(sessionInfo()),paste('mechanism',mechanism,'n',n,
  'replicate',replicate,'seed',seed),paste('PID',Sys.getpid())),file.path(outdir,'session.txt'))
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),
  'Independent local runtime diagnosis. Excluded from all cloud replication counts.'),
  file.path(outdir,'started-at.txt'))
started<-proc.time()[3L];warnings<-character()
profile_run<-function() {
  Rprof(file.path(outdir,'Rprof.out'),interval=.05)
  on.exit(Rprof(NULL))
  run_dataset(mechanism,n,replicate,seed,all_correct=TRUE,dose_max=3L,
    learner_groups=3L,balance_measurement=TRUE,visits='missing_both')
}
cat('Profiling frozen dataset:',mechanism,'n',n,'replication',replicate,
  'seed',seed,'PID',Sys.getpid(),'\n')
result<-tryCatch(withCallingHandlers(profile_run(),warning=function(w) {
  warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
}),error=function(e)list(design=list(mechanism=mechanism,n=n,replicate=replicate,seed=seed),
  errors=data.table(error=conditionMessage(e))))
result$seconds<-proc.time()[3L]-started;result$job_warnings<-unique(warnings)
saveRDS(result,result_file)
profile<-summaryRprof(file.path(outdir,'Rprof.out'))
for(kind in c('by.total','by.self')) {
  table<-data.table(function_name=rownames(profile[[kind]]),profile[[kind]])
  fwrite(table,file.path(outdir,paste0('runtime-',gsub('\\.','-',kind),'.csv')))
}
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),paste('elapsed_seconds',result$seconds),
  paste('final_error_records',nrow(result$errors))),file.path(outdir,'finished-at.txt'))
cat('Runtime diagnosis saved. Elapsed seconds:',result$seconds,
  'final errors:',nrow(result$errors),'\n')
