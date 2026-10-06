# Compare already fitted checkpoints; never combine platforms or refit models.
library(data.table)
cloud<-Sys.getenv('CLOUD_RESULTS');local<-Sys.getenv('LOCAL_RESULTS')
output<-Sys.getenv('COMPARISON_OUTPUT')
stopifnot(nzchar(cloud),nzchar(local),nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
version_file<-Sys.getenv('CLOUD_DEPENDENCY_VERSIONS')
if(nzchar(version_file)) {
  versions<-fread(version_file)
  packages<-c('cmbridge','lmtp','glmnet','SuperLearner','earth','Matrix','data.table','digest')
  differences<-data.table(package=packages,
    local_version=vapply(packages,function(p)as.character(packageVersion(p)),''),
    cloud_version=as.character(package_version(versions$Version[match(packages,versions$Package)])))
  differences[,same_version:=local_version==cloud_version]
  fwrite(differences,file.path(output,'dependency-versions.csv'))
}
rows<-function_rows<-checks<-list()
for(path in list.files(file.path(cloud,'jobs'),'\\.rds$',full.names=TRUE)) {
  other<-file.path(local,'jobs',basename(path))
  if(!file.exists(other))next
  x<-readRDS(path);y<-readRDS(other)
  stopifnot(identical(x$design$seed,y$design$seed),identical(x$folds,y$folds),
    identical(x$learner_folds,y$learner_folds),
    setequal(names(x$fitted_functions),names(y$fitted_functions)))
  checks[[length(checks)+1L]]<-data.table(job=basename(path),seed=x$design$seed,
    outer_splits_equal=identical(x$folds,y$folds),learner_splits_equal=identical(x$learner_folds,y$learner_folds),
    training_combination_counts_equal=isTRUE(all.equal(x$occupancy,y$occupancy)),
    cloud_seconds=x$seconds,local_seconds=y$seconds,
    cloud_fingerprint=x$fingerprint,local_fingerprint=y$fingerprint)
  estimates<-merge(x$results[,.(estimator,horizon,scenario,cloud_estimate=estimate,cloud_se=se)],
    y$results[,.(estimator,horizon,scenario,local_estimate=estimate,local_se=se)],
    by=c('estimator','horizon','scenario'))
  estimates[,`:=`(job=basename(path),seed=x$design$seed,
    estimate_difference=cloud_estimate-local_estimate)]
  rows[[length(rows)+1L]]<-estimates
  for(key in names(x$fitted_functions))for(component in names(x$fitted_functions[[key]])) {
    a<-unlist(x$fitted_functions[[key]][[component]],use.names=FALSE)
    b<-unlist(y$fitted_functions[[key]][[component]],use.names=FALSE)
    stopifnot(length(a)==length(b),all(is.finite(a)),all(is.finite(b)))
    function_rows[[length(function_rows)+1L]]<-data.table(job=basename(path),seed=x$design$seed,
      fit=key,component=component,values=length(a),max_absolute_difference=max(abs(a-b)),
      rmse_difference=sqrt(mean((a-b)^2)))
  }
}
if(length(checks)) {
  fwrite(rbindlist(checks),file.path(output,'matching-checkpoints.csv'))
  fwrite(rbindlist(rows),file.path(output,'estimate-differences.csv'))
  fwrite(rbindlist(function_rows),file.path(output,'fitted-function-differences.csv'))
  print(rbindlist(rows))
} else writeLines('No matching completed seeds yet.',file.path(output,'STATUS.txt'))
writeLines(c('Only saved fits are compared; no models are refitted.',
  'Cloud and Mac results have separate fingerprints and are not pooled.',
  'Equal splits and training counts do not imply numerically identical fitted functions.',
  'Both platforms use R 4.5.2 but different operating systems and BLAS implementations.',
  'The prepared Linux library also has newer glmnet and SuperLearner versions than the local library.',
  'The cmbridge and modified lmtp source archives and installed versions are identical.',
  'Consult each session and dependency-version record before attributing a difference to a particular library.'),
  file.path(output,'NOTES.txt'))
