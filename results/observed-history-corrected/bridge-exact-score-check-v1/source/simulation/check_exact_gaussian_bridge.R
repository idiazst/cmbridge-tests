# Isolated estimated-model diagnostic using the existing exact Gaussian U score.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
args<-commandArgs(TRUE);mechanism<-args[1L];root<-Sys.getenv('EXACT_BRIDGE_OUTPUT')
reference<-Sys.getenv('EXACT_BRIDGE_REFERENCE');stopifnot(nzchar(root),nzchar(reference),
 mechanism%in%c('binary_longitudinal','discrete_dose'))
seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
x<-readRDS(file.path(reference,paste0(mechanism,'-n4000-seed',seed),'result.rds'))
set.seed(seed);g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
require_valid_mechanism(g);d<-draw_data(4000L,g)
prepared<-study_task(d,g,folds=x$folds,learner_folds=x$learner_folds,learner_groups=3L)
task<-prepared$task;stopifnot(identical(d,x$data),identical(task$folds,x$folds),
 identical(task$learner_folds,x$learner_folds))
out<-file.path(root,mechanism);dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(length(list.files(out,'fit-fold[1-3][.]rds$')))stop('Preserve the original saved exact-score fits.')
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),paste('PID',Sys.getpid()),
 'Sieve full conditioning cells; Landweber degree-one conditioning and target; exact Gaussian U scoring.'),file.path(out,'STARTED.txt'))
metrics<-points<-weights<-penalties<-scores<-solvers<-list()
for(fold in 1:3) {
 train<-task$folds[[fold]]$training_set;H<-task$vars$history('A',2L);A<-task$vars$A[[2L]]
 B<-as.matrix(prepared$encoded[train,c(H,A),drop=FALSE])
 V<-as.matrix(cbind(prepared$encoded[train,H,drop=FALSE],d[train,c('C3_covariate','Y3'),drop=FALSE]))
 M<-d$R3[train]
 nested<-lmtp:::make_bridge_nested_folds(task$id[train],d[train,c('R2','R3'),drop=FALSE],task$learner_folds[[fold]])
 library<-cm_library(FALSE,TRUE)
 library$sieve_md$control$instrument_basis<-'cell'
 library$landweber$control$instrument_basis<-'poly'
 library$landweber$control$penalty_scales<-NULL
 for(name in names(library))if(!is.null(library[[name]]$control$penalty_scales)) {
  library[[name]]$control$penalty_ids<-task$id[train]
  library[[name]]$control$penalty_folds<-nested
 }
 warnings<-character();started<-proc.time()[3L]
 model<-withCallingHandlers(cmbridge::fit_bridge_ensemble(B,V,M,library,
   fold_id=task$learner_folds[[fold]],kernel='rbf',kernel_control=list(approximation='exact')),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 elapsed<-proc.time()[3L]-started
 saveRDS(model,file.path(out,paste0('fit-fold',fold,'.rds')))
 writeLines(warnings,file.path(out,paste0('warnings-fold',fold,'.txt')))
 # The score below independently constructs the entire validation kernel.
 gram<-matrix(0,3L,3L)
 for(f in sort(unique(model$fold_id))) {
  valid<-which(model$fold_id==f);training<-which(model$fold_id!=f)
  spec<-cmbridge:::.ensemble_kernel(B[training,,drop=FALSE],'rbf',model$scoring_kernel$control,1L+f)
  stopifnot(is.null(spec$nystrom),spec$control$approximation=='exact')
  z<-sweep(sweep(B[valid,,drop=FALSE],2L,spec$scale$center,'-'),2L,spec$scale$scale,'/')
  squared<-outer(rowSums(z^2),rowSums(z^2),'+')-2*tcrossprod(z)
  K<-exp(-pmax(squared,0)/(2*spec$bandwidth^2));diag(K)<-0
  r<-model$cv_residuals[valid,,drop=FALSE]
  direct<-crossprod(r,K%*%r)/(length(valid)*(length(valid)-1L))
  stopifnot(max(abs(direct-model$fold_gram[[f]]))<1e-10)
  gram<-gram+length(valid)/model$n_scored*direct
 }
 stopifnot(max(abs(gram-model$raw_gram))<1e-10,
   min(eigen(model$gram,symmetric=TRUE)$values)>-1e-10,
   identical(model$fold_id,task$learner_folds[[fold]]),length(model$candidate_failures)==0L)
 cv_list<-c(list(list(sieve_md=model$candidates$sieve_md$penalty_cv)),model$fold_penalty_cv)
 for(scope in seq_along(cv_list))for(name in names(cv_list[[scope]])) {
  cv<-cv_list[[scope]][[name]];if(is.null(cv))next
  chosen<-which(cv$selected)
  stopifnot(length(chosen)==1L,chosen>1L,chosen<nrow(cv),all(cv$scale>0),
   !cv$failed[chosen],is.finite(cv$loss[chosen]),cv$loss[chosen]==min(cv$loss))
  penalties[[length(penalties)+1L]]<-data.table(mechanism,fold,training_scope=scope,
    candidate=name,selected=cv$scale[chosen],loss=cv$loss[chosen],failed_trials=sum(cv$failed),interior=TRUE)
 }
 scores[[length(scores)+1L]]<-data.table(mechanism,fold,elapsed_seconds=elapsed,
  raw_gram_maximum_difference=max(abs(gram-model$raw_gram)),
  psd_adjustment=model$psd_projection$adjustment_norm,candidate_failures=length(model$candidate_failures))
 weights[[length(weights)+1L]]<-data.table(mechanism,fold,candidate=names(model$weights),weight=as.numeric(model$weights))
 for(name in names(model$candidates)) {
  candidate<-model$candidates[[name]]
  solvers[[length(solvers)+1L]]<-data.table(mechanism,fold,candidate=name,
    tolerance_reached=if(name=='landweber')isTRUE(candidate$solver$stopped_by_tolerance) else isTRUE(candidate$solver$converged))
 }
 # Population functions enter only after fitting and all selection.
 p<-x$population;pt<-study_task(p,g)$task;observed<-p$R3==1;included<-p$R2[observed]==1
 Vp<-as.matrix(cbind(pt$natural[observed,H,drop=FALSE],p[observed,c('C3_covariate','Y3'),drop=FALSE]))
 truth<-population_truth_functions(p,g)$beta[observed,2L]
 keys<-cell_key(pt$natural[,c(H,A),drop=FALSE])
 for(name in names(c(list(ensemble=model),model$candidates))) {
  candidate<-if(name=='ensemble')model else model$candidates[[name]]
  fitted<-as.numeric(predict(candidate,Vp));stopifnot(all(is.finite(fitted)),all(fitted>=1))
  residual<-rep(-1,nrow(p));residual[observed]<-fitted-1
  equation<-data.table(cell=keys,probability=p$probability,residual=residual)
  equation<-equation[,.(probability=sum(probability),error=sum(probability*residual)/sum(probability)),by=cell]
  metrics[[length(metrics)+1L]]<-data.table(mechanism,fold,candidate=name,
    bridge_rmse=sqrt(sum(p$probability[observed][included]*(fitted[included]-truth[included])^2)/sum(p$probability[observed][included])),
    equation_rmse=sqrt(sum(equation$probability*equation$error^2)),
    equation_maximum_error=max(abs(equation$error)),minimum=min(fitted),maximum=max(fitted),
    application_bound_n=sum(fitted<=-100|fitted>=100))
  points[[length(points)+1L]]<-data.table(mechanism,fold,candidate=name,
    reference=truth[included],fitted=fitted[included],probability=p$probability[observed][included])
 }
 for(pair in list(list(metrics,'bridge-errors.csv'),list(points,'plot-points.csv'),
    list(weights,'ensemble-weights.csv'),list(penalties,'penalty-audit.csv'),list(scores,'score-audit.csv'),list(solvers,'solver-audit.csv')))
  fwrite(rbindlist(pair[[1L]]),file.path(out,pair[[2L]]))
 cat(format(Sys.time()),'Saved exact-score outer training sample',fold,'seconds',elapsed,'\n')
}
writeLines(c(format(Sys.time(),tz='UTC',usetz=TRUE),'All three bridge-only checks complete; complete-estimator performance not assessed.'),file.path(out,'FINISHED.txt'))
