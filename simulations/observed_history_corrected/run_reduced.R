study_source <- Sys.getenv('STUDY_SOURCE', 'simulations/observed_history_corrected')
source(file.path(study_source, 'study.R'))
outdir <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/missing-both-study-v2')
dir.create(file.path(outdir, 'jobs'), recursive = TRUE, showWarnings = FALSE)
workers <- as.integer(Sys.getenv('SIM_WORKERS', '8'))
replications <- as.integer(Sys.getenv('SIM_REPS', '200'))
sizes <- as.integer(strsplit(Sys.getenv('SIM_SIZES', '500,1000,4000'), ',')[[1]])
mechanisms <- c('binary_longitudinal', 'discrete_dose')
stopifnot(workers >= 1L, replications >= 1L, all(sizes %in% c(500L,1000L,4000L)))
grid <- CJ(mechanism = mechanisms, n = sizes, replicate = seq_len(replications))
grid[, seed := 5000000L + match(mechanism, mechanisms) * 100000L +
       match(n, c(500L,1000L,4000L)) * 1000L + replicate]
grid[, job := sprintf('%s-n%d-r%03d', mechanism, n, replicate)]
lmtp_source <- Sys.getenv('LMTP_SOURCE', '../lmtp')
cmbridge_source <- Sys.getenv('CMBRIDGE_SOURCE', '../cmbridge')
sources <- c(file.path(study_source, c('dgp.R','study.R','run_reduced.R')),
  list.files(file.path(lmtp_source,'R'),'\\.R$',full.names=TRUE),
  list.files(file.path(cmbridge_source,'R'),'\\.R$',full.names=TRUE),
  file.path(lmtp_source,'DESCRIPTION'),file.path(cmbridge_source,'DESCRIPTION'))
fingerprint <- digest::digest(c(tools::md5sum(sort(sources)),packageVersion('lmtp'),packageVersion('cmbridge')))
manifest <- file.path(outdir,'source-fingerprint.txt')
if (file.exists(manifest) && readLines(manifest)[1] != fingerprint)
  stop('Source differs from this run. Preserve its checkpoints and use a new output directory.')
writeLines(fingerprint,manifest)
started_file<-file.path(outdir,'run-started-at.txt')
if(!file.exists(started_file))writeLines(format(Sys.time(),tz='America/New_York',usetz=TRUE),started_file)
fwrite(grid,file.path(outdir,'design.csv'))
completed_initial <- length(list.files(file.path(outdir,'jobs'),'\\.rds$'))
writeLines(c('Running the revised study with missingness at both follow-ups.',
  paste('Completed datasets:',completed_initial,'of',nrow(grid)),
  paste('Remaining datasets:',nrow(grid)-completed_initial),paste('Workers:',workers)),
  file.path(outdir,'STATUS.txt'))
writeLines(c(paste('Completed datasets:',completed_initial,'of',nrow(grid)),
  paste('Remaining datasets:',nrow(grid)-completed_initial),paste('Workers:',workers),
  paste('Updated:',format(Sys.time(),tz='America/New_York',usetz=TRUE))),
  file.path(outdir,'PROGRESS.txt'))
fwrite(data.table(estimator=c('sdr','tmle'),beta_correct=TRUE,lambda_correct=TRUE,
  ratio_correct=TRUE,m_correct=TRUE),file.path(outdir,'specifications.csv'))
saveRDS(list(workers=workers,replications=replications,sizes=sizes,mechanisms=mechanisms,
  dose_levels=0:3,all_correct=TRUE,learner_groups=3L,measurement_balanced=TRUE,
  visits='missing_both'),file.path(outdir,'design-details.rds'))
writeLines(c(capture.output(sessionInfo()),paste('workers',workers),paste('fingerprint',fingerprint)),
  file.path(outdir,'session.txt'))
selected <- order(grid$replicate,grid$n,grid$mechanism)
one <- function(i) {
  config<-grid[i]; path<-file.path(outdir,'jobs',paste0(config$job,'.rds'))
  if (file.exists(path)) {
    previous<-readRDS(path)
    stopifnot(identical(previous$fingerprint,fingerprint),previous$design$seed==config$seed)
    return(invisible(TRUE))
  }
  cat('Starting',config$job,'worker',Sys.getpid(),'\n')
  started<-proc.time()[3]; warnings<-character()
  result<-tryCatch(withCallingHandlers(
    run_dataset(config$mechanism,config$n,config$replicate,config$seed,all_correct=TRUE,dose_max=3L,learner_groups=3L,balance_measurement=TRUE,visits='missing_both'),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e)list(design=as.list(config[,.(mechanism,n,replicate,seed)]),
      errors=data.table(error=conditionMessage(e))))
  result$fingerprint<-fingerprint; result$seconds<-proc.time()[3]-started
  result$job_warnings<-unique(warnings)
  tmp<-paste0(path,'.tmp');saveRDS(result,tmp,compress='gzip');stopifnot(file.rename(tmp,path))
  status<-if(nrow(result$errors))paste(unique(result$errors$error),collapse=' | ')else 'OK'
  cat(config$job,sprintf('%.1fs',result$seconds),status,'\n');invisible(TRUE)
}
# Several replications per batch keep all eight workers occupied across six design cells.
batches<-split(selected,ceiling(seq_along(selected)/(workers*3L)))
for (b in seq_along(batches)) {
  batch<-batches[[b]]
  answers<-parallel::mclapply(batch,one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  saved<-file.exists(file.path(outdir,'jobs',paste0(grid$job[batch],'.rds')))
  if(length(answers)!=length(batch)||!all(vapply(answers,isTRUE,FALSE))||!all(saved))
    stop('A worker did not save its checkpoint; inspect and resume unfinished jobs.')
  completed_paths<-list.files(file.path(outdir,'jobs'),'\\.rds$')
  writeLines(c(paste('Completed datasets:',length(completed_paths),'of',nrow(grid)),
    paste('Remaining datasets:',nrow(grid)-length(completed_paths)),
    paste('Workers:',workers),paste('Updated:',format(Sys.time(),tz='America/New_York',usetz=TRUE)),
    'Coverage and convergence conclusions remain preliminary until all planned replications finish.'),file.path(outdir,'PROGRESS.txt'))
  if(Sys.getenv('SIM_POSTPROCESS','0')=='1' && (b==1L || b%%5L==0L || b==length(batches))) {
    for(script in c('summarize.R','write_reduced_report.R')) {
      status<-system2(file.path(R.home('bin'),'Rscript'),shQuote(file.path(study_source,script)))
      if(status!=0)stop('Postprocessing failed; saved simulation checkpoints remain available.')
    }
  }
}
writeLines(format(Sys.time(),tz='America/New_York',usetz=TRUE),file.path(outdir,'run-finished-at.txt'))
cat('Requested reduced-study datasets completed.\n')
