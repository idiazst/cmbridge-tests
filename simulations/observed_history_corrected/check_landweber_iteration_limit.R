# Same estimated-model bridge, same positive CV choice; solver sensitivity only.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
root<-Sys.getenv('LANDWEBER_ITERATION_ROOT');base<-Sys.getenv('LANDWEBER_ITERATION_BASE')
mechanism<-Sys.getenv('LANDWEBER_ITERATION_MECHANISM','binary_longitudinal')
fold<-as.integer(Sys.getenv('LANDWEBER_ITERATION_FOLD','1'))
seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
folder<-paste0(mechanism,'-n4000-seed',seed)
stopifnot(nzchar(root),nzchar(base));dir.create(root,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(root,'fit.rds')))stop('Preserve the original iteration check.')
x<-readRDS(file.path(base,'ensemble-stability-v6',folder,'result.rds'))
fit_file<-Sys.getenv('LANDWEBER_ITERATION_FIT')
m<-if(nzchar(fit_file))readRDS(fit_file)$beta_fits[[2L]]$candidates$landweber else
 readRDS(file.path(base,'landweber-weight-cv-v1/binary_longitudinal/bridge-fit.rds'))$candidates$landweber
stopifnot(!isTRUE(m$solver$stopped_by_tolerance),m$tuning$weight_ridge>0)
g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
task<-study_task(x$data,g,folds=x$folds,learner_folds=x$learner_folds,learner_groups=3L)$task
rows<-x$folds[[fold]]$training_set;H<-task$vars$history('A',2L);A<-task$vars$A[[2L]]
B<-as.matrix(task$natural[rows,c(H,A),drop=FALSE])
V<-as.matrix(cbind(task$natural[rows,H,drop=FALSE],x$data[rows,c('C3_covariate','Y3'),drop=FALSE]))
control<-m$tuning
for(name in c('step','maximum_step','iterations_used','final_delta'))control[[name]]<-NULL
control$max_iter<-10000L
started<-proc.time()[3L];fit<-cmbridge::fit_bridge(B,V,x$data$R3[rows],'landweber',control)
elapsed<-proc.time()[3L]-started;saveRDS(fit,file.path(root,'fit.rds'))
p<-x$population;pt<-study_task(p,g)$task;observed<-p$R3==1
Vp<-as.matrix(cbind(pt$natural[observed,H,drop=FALSE],p[observed,c('C3_covariate','Y3'),drop=FALSE]))
old<-predict(m,Vp);new<-predict(fit,Vp)
truth<-population_truth_functions(p,g)$beta[observed,2L]
probability<-p$probability[observed];included<-p$R2[observed]==1
rmse<-function(value)sqrt(sum(probability[included]*(value[included]-truth[included])^2)/sum(probability[included]))
fwrite(data.table(mechanism=mechanism,fold=fold,previous_bridge_rmse=rmse(old),new_bridge_rmse=rmse(new),maximum_prediction_change=max(abs(new-old)),
 elapsed_seconds=elapsed,previous_iterations=m$solver$iterations,
 new_iterations=fit$solver$iterations,previous_gradient=m$solver$coefficient_gradient,
 new_gradient=fit$solver$coefficient_gradient,previous_tolerance_reached=m$solver$stopped_by_tolerance,
 new_tolerance_reached=fit$solver$stopped_by_tolerance,weight_ridge=fit$tuning$weight_ridge,
 tol=fit$tuning$tol),file.path(root,'comparison.csv'))
print(fread(file.path(root,'comparison.csv')))
