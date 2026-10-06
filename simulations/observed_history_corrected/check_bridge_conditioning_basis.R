# A separate estimated-model diagnostic on one saved n=4000 training sample.
# The frozen full study is never edited or pooled with these alternative fits.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
input<-Sys.getenv('BRIDGE_CONDITIONING_INPUT');output<-Sys.getenv('BRIDGE_CONDITIONING_OUTPUT')
basis<-Sys.getenv('BRIDGE_CONDITIONING_BASIS')
stopifnot(nzchar(input),nzchar(output),basis%in%c('polynomial','joint_categories'))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(output,'bridge-fit.rds')))stop('Preserve the saved fit; do not overwrite it.')
x<-readRDS(input);config<-x$design
stopifnot(config$n==4000L,config$mechanism=='binary_longitudinal',config$seed==5103002L)
set.seed(config$seed);g<-make_mechanism(config$mechanism,dose_max=3L,visits='missing_both')
require_valid_mechanism(g);d<-draw_data(config$n,g)
prepared<-study_task(d,g,folds=3L,learner_groups=3L,balance_measurement=TRUE)
task<-prepared$task
stopifnot(identical(task$folds,x$folds),identical(task$learner_folds,x$learner_folds))
fold<-1L;s<-2L;train<-task$folds[[fold]]$training_set
H<-task$vars$history('A',s);A<-task$vars$A[[s]]
B<-as.matrix(prepared$encoded[train,c(H,A),drop=FALSE])
V<-as.matrix(cbind(prepared$encoded[train,H,drop=FALSE],
  d[train,c('C3_covariate','Y3'),drop=FALSE]))
M<-d$R3[train]
nested<-lmtp:::make_bridge_nested_folds(task$id[train],d[train,c('R2','R3'),drop=FALSE],
  task$learner_folds[[fold]])
library<-cm_library(FALSE,TRUE)
if(basis=='joint_categories')for(name in c('sieve_md','landweber'))
  library[[name]]$control$instrument_basis<-'cell'
for(name in names(library))if(!is.null(library[[name]]$control$penalty_scales)) {
  library[[name]]$control$penalty_ids<-task$id[train]
  library[[name]]$control$penalty_folds<-nested
}
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),paste('PID',Sys.getpid()),
  paste('basis',basis,'seed',config$seed,'outer_training_sample',fold),
  'Separate diagnostic; no known-function fits or zero-penalty fits.'),file.path(output,'STARTED.txt'))
warnings<-character();started<-proc.time()[3L]
model<-tryCatch(withCallingHandlers(cmbridge::fit_bridge_ensemble(B,V,M,
  library=library,fold_id=task$learner_folds[[fold]],kernel='rbf'),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e){writeLines(conditionMessage(e),file.path(output,'ERROR.txt'));stop(e)})
elapsed<-proc.time()[3L]-started
saveRDS(model,file.path(output,'bridge-fit.rds'))
writeLines(warnings,file.path(output,'warnings.txt'))
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),paste('elapsed_seconds',elapsed)),
  file.path(output,'FINISHED.txt'))
p<-x$population;pt<-study_task(p,g)$task
measured<-p$R3==1
Vp<-as.matrix(cbind(pt$natural[measured,H,drop=FALSE],
  p[measured,c('C3_covariate','Y3'),drop=FALSE]))
truth<-population_truth_functions(p,g)$beta[measured,s]
predictions<-list(ensemble=as.numeric(predict(model,Vp)))
for(name in names(model$candidates))if(!is.null(model$candidates[[name]]))
  predictions[[name]]<-as.numeric(predict(model$candidates[[name]],Vp))
baseline<-x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',fold,sep='/')]]$beta[measured,s]
if(basis=='polynomial')stopifnot(max(abs(predictions$ensemble-baseline))<1e-9)
keys<-cell_key(pt$natural[,c(H,A),drop=FALSE]);records<-list();points<-list()
for(name in names(predictions)) {
  fitted<-predictions[[name]];residual<-rep(-1,nrow(p));residual[measured]<-fitted-1
  cells<-data.table(cell=keys,probability=p$probability,residual=residual)
  cells<-cells[,.(probability=sum(probability),error=sum(probability*residual)/sum(probability)),by=cell]
  selected<-p$R2[measured]==1
  records[[length(records)+1L]]<-data.table(basis=basis,candidate=name,
    bridge_rmse_R2_1=sqrt(sum(p$probability[measured][selected]*(fitted[selected]-truth[selected])^2)/
      sum(p$probability[measured][selected])),
    defining_equation_rmse=sqrt(sum(cells$probability*cells$error^2)),
    defining_equation_max_error=max(abs(cells$error)),
    raw_minimum=min(fitted),raw_maximum=max(fitted),
    application_bound_fraction=mean(fitted<=-100|fitted>=100))
  points[[length(points)+1L]]<-data.table(basis=basis,candidate=name,
    probability=p$probability[measured][selected],reference=truth[selected],fitted=fitted[selected])
}
fwrite(rbindlist(records),file.path(output,'population-bridge-checks.csv'))
fwrite(rbindlist(points),file.path(output,'bridge-plot-points.csv'))
fwrite(data.table(candidate=names(model$weights),weight=as.numeric(model$weights)),
  file.path(output,'ensemble-weights.csv'))
for(name in names(model$candidates)) {
  candidate<-model$candidates[[name]]
  if(!is.null(candidate$penalty_cv))fwrite(candidate$penalty_cv,
    file.path(output,paste0(name,'-penalty-cv.csv')))
}
writeLines(c(paste('elapsed_seconds',elapsed),paste('candidate_failure_records',length(model$candidate_failures)),
  paste('baseline_prediction_max_difference',max(abs(predictions$ensemble-baseline)))),
  file.path(output,'AUDIT.txt'))
print(rbindlist(records));print(model$weights)
