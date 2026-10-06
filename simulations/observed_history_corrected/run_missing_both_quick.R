# One n=10000 dataset per mechanism. This does not launch the repeated study.
study_source <- Sys.getenv('STUDY_SOURCE', 'simulations/observed_history_corrected')
source(file.path(study_source, 'study.R'))
out <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/missing-both-n10000')
dir.create(file.path(out,'jobs'),recursive=TRUE,showWarnings=FALSE)
design <- data.table(mechanism=c('binary_longitudinal','discrete_dose'),
  n=10000L,replicate=1L,seed=c(6100001L,6200001L))
fwrite(design,file.path(out,'design.csv'))
writeLines(c(capture.output(sessionInfo()),paste('library',.libPaths()[1L])),file.path(out,'session.txt'))
writeLines(format(Sys.time(),tz='America/New_York',usetz=TRUE),file.path(out,'run-started-at.txt'))
one <- function(i) {
  config <- design[i]; path <- file.path(out,'jobs',paste0(config$mechanism,'.rds'))
  if(file.exists(path))stop('Existing checkpoint: preserve it and use a new output directory.')
  cat('Starting',config$mechanism,'n=10000 seed=',config$seed,'\n')
  started <- proc.time()[3L]; warnings <- character()
  result <- tryCatch(withCallingHandlers(
    run_dataset(config$mechanism,config$n,config$replicate,config$seed,
      all_correct=TRUE,dose_max=3L,learner_groups=3L,balance_measurement=TRUE,
      visits='missing_both',store_scores=TRUE),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e)list(design=as.list(config),errors=data.table(error=conditionMessage(e)),
      call=deparse(conditionCall(e))))
  result$seconds <- proc.time()[3L]-started
  result$job_warnings <- unique(warnings)
  tmp <- paste0(path,'.tmp');saveRDS(result,tmp,compress='gzip')
  stopifnot(file.rename(tmp,path))
  cat('Finished',config$mechanism,sprintf('%.1fs',result$seconds),
    if(nrow(result$errors))paste(result$errors$error,collapse=' | ')else 'OK','\n')
  TRUE
}
answers <- parallel::mclapply(seq_len(nrow(design)),one,mc.cores=2L,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(all(vapply(answers,isTRUE,FALSE)))
jobs <- lapply(file.path(out,'jobs',paste0(design$mechanism,'.rds')),readRDS)
errors <- rbindlist(lapply(seq_along(jobs),function(i){
  e<-copy(jobs[[i]]$errors);if(nrow(e))e[,mechanism:=design$mechanism[i]];e}),fill=TRUE)
if(nrow(errors))fwrite(errors,file.path(out,'errors.csv'))
results <- rbindlist(lapply(jobs,`[[`,'results'),fill=TRUE)
components <- rbindlist(lapply(jobs,`[[`,'eif_components'),fill=TRUE)
if(nrow(results))fwrite(results,file.path(out,'estimates.csv'))
if(nrow(components)) {
  fwrite(components,file.path(out,'person-eif-contributions.csv.gz'))
  summary <- components[,.(
    sequential_estimate=mean(sequential_score),
    bridge_adjoint_correction=mean(bridge_adjoint),
    final_estimate=mean(sequential_score+bridge_adjoint),
    sequential_eif_mean=mean(sequential_eif),
    bridge_adjoint_eif_mean=mean(bridge_adjoint),total_eif_mean=mean(total_eif),
    sequential_eif_sd=sd(sequential_eif),bridge_adjoint_eif_sd=sd(bridge_adjoint),
    eif_covariance=cov(sequential_eif,bridge_adjoint),total_eif_sd=sd(total_eif),
    sequential_component_se=sd(sequential_eif)/sqrt(.N),
    bridge_adjoint_component_se=sd(bridge_adjoint)/sqrt(.N),
    total_se=sd(total_eif)/sqrt(.N),persons=.N),by=.(mechanism,estimator,horizon)]
  stopifnot(max(abs(summary$total_eif_mean))<1e-12,
    max(abs(summary$total_eif_sd^2-summary$sequential_eif_sd^2-
      summary$bridge_adjoint_eif_sd^2-2*summary$eif_covariance))<1e-10)
  fwrite(summary,file.path(out,'eif-summary.csv'))
}
diagnostics <- rbindlist(lapply(jobs,`[[`,'diagnostics'),fill=TRUE)
if(nrow(diagnostics)) {
  fwrite(diagnostics,file.path(out,'population-diagnostics-by-fold.csv'))
  summary <- diagnostics[,.(population_bias=sum(fold_weight*population_bias),
    sequential_remainder=sum(fold_weight*sequential_remainder),
    bridge_remainder=sum(fold_weight*bridge_remainder),
    maximum_bridge_identity_error=max(abs(bridge_identity_error))),
    by=.(mechanism,estimator,horizon)]
  stopifnot(max(abs(summary$population_bias-summary$sequential_remainder-
    summary$bridge_remainder))<1e-10)
  fwrite(summary,file.path(out,'population-remainders.csv'))
}
for(field in c('tuning','ensemble_weights','occupancy')) {
  table <- rbindlist(lapply(seq_along(jobs),function(i){
    v<-copy(jobs[[i]][[field]]);if(nrow(v))v[,mechanism:=design$mechanism[i]];v}),fill=TRUE)
  if(nrow(table))fwrite(table,file.path(out,paste0(field,'.csv')))
}
missing <- rbindlist(lapply(seq_along(jobs),function(i){
  d<-jobs[[i]]$data;if(is.null(d))return(NULL)
  fwrite(d,file.path(out,paste0(design$mechanism[i],'-data.csv.gz')))
  data.table(mechanism=design$mechanism[i],time=1:3,
    sample_missing_percent=100*vapply(d[,paste0('R',1:3)],function(x)mean(x==0),0.0))
}))
if(nrow(missing))fwrite(missing,file.path(out,'sample-missingness.csv'))
writeLines(format(Sys.time(),tz='America/New_York',usetz=TRUE),file.path(out,'run-finished-at.txt'))
if(nrow(errors))stop('A requested estimator failed; its data and error are preserved.')
stopifnot(nrow(results)==8L,nrow(components)==80000L)
writeLines('Both single datasets and both estimators completed. The repeated-sample study remains stopped.',file.path(out,'STATUS.txt'))
cat('Both requested single-dataset checks completed.\n')
