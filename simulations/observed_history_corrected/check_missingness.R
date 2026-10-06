# Diagnostic only: reproduce generated samples from saved seeds; do not fit.
library(data.table)
outdir <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/ensemble-study')
study_source <- file.path(outdir, 'source', 'simulations', 'observed_history_corrected')
source(file.path(study_source, 'dgp.R'))
progress <- readLines(file.path(outdir, 'PROGRESS.txt'))
completed <- as.integer(sub('Completed datasets: ([0-9]+).*', '\\1', progress[1L]))
max_rep <- completed %/% 6L
design <- fread(file.path(outdir, 'design.csv'))[replicate <= max_rep]
paths <- file.path(outdir, 'jobs', paste0(design$job, '.rds'))
stopifnot(all(file.exists(paths)))
check <- file.path(outdir, 'missingness-check')
dir.create(check, showWarnings=FALSE)
mechanisms <- setNames(lapply(unique(design$mechanism), function(name)
  make_mechanism(name, dose_max=3L)), unique(design$mechanism))
expected <- rbindlist(lapply(mechanisms, function(g) data.table(
  mechanism=g$mechanism, time=1:3,
  missing_percent=100*as.numeric(apply(g$natural_states, 2L, function(x) sum(x[,1L]))),
  policy_missing_percent=100*as.numeric(apply(g$policy_states, 2L, function(x) sum(x[,1L]))))))
conditional <- list()
for (g in mechanisms) {
  cells <- list()
  for (b in 1:2) for(h in seq_len(nrow(g$L))) for(a in seq_len(nrow(g$A))) {
    w <- g$natural_states[b,2L,h]*g$g[b,2L,h,a]*g$health[b,2L,h,a,]
    cells[[length(cells)+1L]] <- data.table(outcome=g$C$outcome,
      probability=w, missing_probability=w*(1-g$measurement[b,2L,h,]))
  }
  values <- rbindlist(cells)[,.(missing_percent=100*sum(missing_probability)/sum(probability)),by=outcome]
  values[,mechanism:=g$mechanism]
  conditional[[length(conditional)+1L]] <- values
}
sample_rows <- diagnoses <- estimates <- occupancy <- list()
for (i in seq_len(nrow(design))) {
  config <- design[i]
  g <- mechanisms[[config$mechanism]]
  set.seed(config$seed)
  d <- draw_data(config$n,g)
  for(t in 1:3) {
    miss <- d[[paste0('R',t)]]==0
    stopifnot(identical(miss,is.na(d[[paste0('Y',t)]])),
              identical(miss,is.na(d[[paste0('C',t,'_covariate')]])))
    sample_rows[[length(sample_rows)+1L]] <- data.table(mechanism=config$mechanism,
      n=config$n,replicate=config$replicate,time=t,missing=sum(miss),people=nrow(d))
  }
  x <- readRDS(paths[i])
  if(nrow(x$diagnostics))diagnoses[[length(diagnoses)+1L]] <- x$diagnostics
  if(nrow(x$results))estimates[[length(estimates)+1L]] <- x$results
  if(nrow(x$occupancy)) {
    oc <- copy(x$occupancy[kind=='V' & measured>0])
    oc <- oc[,.(represented_measured_cells=.N,cells_under_5=sum(measured<5),
                cells_under_10=sum(measured<10),measured_training_people=sum(measured)),
             by=.(fold,horizon)]
    oc[,`:=`(mechanism=config$mechanism,n=config$n,replicate=config$replicate)]
    occupancy[[length(occupancy)+1L]] <- oc
  }
}
ss <- rbindlist(sample_rows)
empirical <- ss[,.(datasets=.N,missing=sum(missing),people=sum(people),
                  missing_percent=100*sum(missing)/sum(people)),by=.(mechanism,time)]
by_size <- ss[,.(datasets=.N,missing_percent=100*sum(missing)/sum(people)),by=.(mechanism,n,time)]
dg <- rbindlist(diagnoses)
terms <- c('population_bias','sequential_remainder','bridge_remainder')
per_rep <- dg[,lapply(.SD,function(z)sum(z*fold_weight)),
              by=.(mechanism,n,replicate,estimator,horizon),.SDcols=terms]
contributions <- per_rep[,c(list(replications=.N),lapply(.SD,mean)),
                       by=.(mechanism,n,estimator,horizon),.SDcols=terms]
res <- rbindlist(estimates)
performance <- res[,.(replications=.N,truth=mean(truth),bias=mean(estimate-truth),
                rmse=sqrt(mean((estimate-truth)^2)),mean_se=mean(se),
                coverage=mean(lower<=truth & upper>=truth)),by=.(mechanism,n,estimator,horizon)]
oc <- rbindlist(occupancy)[,.(training_fits=.N,
  mean_represented_measured_cells=mean(represented_measured_cells),
  mean_cells_under_5=mean(cells_under_5),mean_cells_under_10=mean(cells_under_10),
  mean_measured_training_people=mean(measured_training_people)),by=.(mechanism,n,horizon)]
fwrite(expected,file.path(check,'expected-missingness.csv'))
fwrite(empirical,file.path(check,'sample-missingness.csv'))
fwrite(by_size,file.path(check,'sample-missingness-by-size.csv'))
fwrite(rbindlist(conditional),file.path(check,'missingness-by-outcome.csv'))
fwrite(contributions,file.path(check,'error-contributions.csv'))
fwrite(performance,file.path(check,'performance.csv'))
fwrite(oc,file.path(check,'measured-cell-counts.csv'))
file.copy('simulations/observed_history_corrected/check_missingness.R',file.path(check,'check_missingness.R'),overwrite=TRUE)
writeLines(c(paste('Completed replication batches:',max_rep),
  paste('Datasets:',nrow(design)),format(Sys.time(),tz='America/New_York',usetz=TRUE),
  'Samples regenerated exactly from the saved seeds and frozen generating mechanisms; no estimators fitted.'),
  file.path(check,'SNAPSHOT.txt'))
print(expected);print(empirical);print(rbindlist(conditional))
print(performance[n==4000L]);print(contributions[n==4000L]);print(oc[n==4000L & horizon==2L])
