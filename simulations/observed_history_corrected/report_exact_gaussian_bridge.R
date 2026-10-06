# Read-only figures and score diagnostics; no model refits or study adoption.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
root<-Sys.getenv('EXACT_BRIDGE_OUTPUT');reference<-Sys.getenv('EXACT_BRIDGE_REFERENCE')
stopifnot(nzchar(root),nzchar(reference))
all_metrics<-all_points<-all_weights<-score_comparisons<-all_cv_paths<-list()
for(mechanism in c('binary_longitudinal','discrete_dose')) {
 out<-file.path(root,mechanism);stopifnot(file.exists(file.path(out,'FINISHED.txt')))
 metrics<-fread(file.path(out,'bridge-errors.csv'));points<-fread(file.path(out,'plot-points.csv'))
 all_metrics[[mechanism]]<-metrics;all_points[[mechanism]]<-points
 all_weights[[mechanism]]<-fread(file.path(out,'ensemble-weights.csv'))
 seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
 x<-readRDS(file.path(reference,paste0(mechanism,'-n4000-seed',seed),'result.rds'))
 g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
 prepared<-study_task(x$data,g,folds=x$folds,learner_folds=x$learner_folds,learner_groups=3L)
 task<-prepared$task
 for(fold in 1:3) {
  model<-readRDS(file.path(out,paste0('fit-fold',fold,'.rds')))
  train<-task$folds[[fold]]$training_set;H<-task$vars$history('A',2L);A<-task$vars$A[[2L]]
  B<-as.matrix(prepared$encoded[train,c(H,A),drop=FALSE])
  nested<-lmtp:::make_bridge_nested_folds(task$id[train],x$data[train,c('R2','R3')],task$learner_folds[[fold]])
  sieve<-model$candidates$sieve_md
  stopifnot(identical(sieve$penalty_ids,task$id[train]),
    identical(sieve$penalty_fold_id,nested(task$id[train])),
    length(intersect(sieve$penalty_ids,task$id[task$folds[[fold]]$validation_set]))==0L)
  cv_scopes<-c(list(list(sieve_md=sieve$penalty_cv)),model$fold_penalty_cv)
  for(scope in seq_along(cv_scopes))for(name in names(cv_scopes[[scope]])) {
   cv<-cv_scopes[[scope]][[name]];if(is.null(cv))next
   path<-as.data.table(copy(cv))
   path[,`:=`(mechanism=mechanism,outer_fold=fold,training_scope=scope,candidate=name)]
   all_cv_paths[[length(all_cv_paths)+1L]]<-path
  }
  approximation_gram<-matrix(0,3L,3L)
  for(f in sort(unique(model$fold_id))) {
   validation<-which(model$fold_id==f);training<-which(model$fold_id!=f)
   control<-model$scoring_kernel$control;control$approximation<-'nystrom'
   spec<-cmbridge:::.ensemble_kernel(B[training,,drop=FALSE],'rbf',control,1L+f)
   z<-sweep(sweep(B[validation,,drop=FALSE],2L,spec$scale$center,'-'),2L,spec$scale$scale,'/')
   centers<-spec$nystrom$centers
   squared<-outer(rowSums(z^2),rowSums(centers^2),'+')-2*tcrossprod(z,centers)
   features<-exp(-pmax(squared,0)/(2*spec$bandwidth^2))%*%spec$nystrom$transform
   K<-tcrossprod(features);diag(K)<-0
   residual<-model$cv_residuals[validation,,drop=FALSE]
   direct<-crossprod(residual,K%*%residual)/(length(validation)*(length(validation)-1L))
   package<-cmbridge:::.ensemble_gram(residual,B[validation,,drop=FALSE],spec)
   stopifnot(max(abs(direct-package))<1e-10)
   approximation_gram<-approximation_gram+length(validation)/model$n_scored*direct
  }
  projected<-cmbridge:::.project_gram_psd(approximation_gram)$gram
  weights<-cmbridge:::.simplex_qp(projected)$weights
  score_comparisons[[length(score_comparisons)+1L]]<-data.table(mechanism,fold,
   maximum_raw_gram_change=max(abs(approximation_gram-model$raw_gram)),
   maximum_weight_change=max(abs(weights-model$weights)),
   exact_sieve_weight=model$weights['sieve_md'],approximate_sieve_weight=weights[1L],
   exact_landweber_weight=model$weights['landweber'],approximate_landweber_weight=weights[2L],
   exact_pmmr_weight=model$weights['pmmr'],approximate_pmmr_weight=weights[3L])
 }
 points[,training_sample:=factor(fold)]
 labels<-c(ensemble='Ensemble',sieve_md='Sieve minimum distance',landweber='Landweber',pmmr='PMMR')
 points[,candidate_label:=factor(labels[candidate],levels=labels)]
 for(ensemble_only in c(FALSE,TRUE)) {
  selected<-if(ensemble_only)points[candidate=='ensemble'] else points
  limits<-range(selected$reference,selected$fitted)
  plot<-ggplot(selected,aes(reference,fitted,colour=training_sample))+
   geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey35')+
   geom_point(alpha=.55,size=1.3)+coord_equal(xlim=limits,ylim=limits)+
   scale_colour_manual(values=c('#2166ac','#b35806','#5e3c99'))+
   labs(x='Population bridge solution',y='Estimated bridge',colour='Training sample',
    title=paste(if(mechanism=='binary_longitudinal')'Binary treatment' else 'Numerical dose',
       'n = 4,000: exact Gaussian U selection'),
    subtitle='Outcome at time 3; intermediate and final visits observed',
    caption='Every combination and all three training samples; equal axes, no clipping.')+
   theme_bw(base_size=11)+theme(legend.position='bottom')
  if(!ensemble_only)plot<-plot+facet_wrap(~candidate_label,ncol=2)
  stem<-paste0(mechanism,if(ensemble_only)'-ensemble' else '-all-candidates')
  for(ext in c('png','pdf'))ggsave(file.path(root,paste0(stem,'.',ext)),plot,
    width=if(ensemble_only)7 else 10,height=if(ensemble_only)6.5 else 9,dpi=180)
 }
}
fwrite(rbindlist(all_metrics),file.path(root,'all-bridge-errors.csv'))
fwrite(rbindlist(all_points),file.path(root,'all-plot-points.csv'))
fwrite(rbindlist(all_weights),file.path(root,'all-ensemble-weights.csv'))
fwrite(rbindlist(score_comparisons),file.path(root,'fixed-candidate-score-comparison.csv'))
fwrite(rbindlist(all_cv_paths,fill=TRUE),file.path(root,'all-penalty-cv-paths.csv'))
summary<-rbindlist(all_metrics)[,.(bridge_rmse=sqrt(mean(bridge_rmse^2)),
 equation_rmse=sqrt(mean(equation_rmse^2))),by=.(mechanism,candidate)]
fwrite(summary,file.path(root,'summary.csv'))
cat('Verified training-only penalty IDs and direct approximate U scores; saved actual bridge-only plots. No refits or complete-estimator claim.\n')
