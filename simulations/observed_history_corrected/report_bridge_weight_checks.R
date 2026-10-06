# Read-only audit of saved estimated-model diagnostics; no nuisance refits.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
base<-Sys.getenv('WEIGHT_CHECK_BASE');stopifnot(nzchar(base))
cv_root<-file.path(base,'landweber-weight-cv-v1')
weight_root<-file.path(base,'bridge-conditioning-weight-check-v1')
metrics<-weights<-points<-audits<-selections<-curves<-solvers<-list()
for(mechanism in c('binary_longitudinal','discrete_dose')) {
 seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
 folder<-paste0(mechanism,'-n4000-seed',seed)
 x<-readRDS(file.path(base,'ensemble-stability-v6',folder,'result.rds'))
 set.seed(seed);g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
 d<-draw_data(4000L,g)
 prepared<-study_task(d,g,folds=x$folds,learner_folds=x$learner_folds,learner_groups=3L)
 stopifnot(identical(d,x$data),identical(prepared$task$folds,x$folds),
   identical(prepared$task$learner_folds,x$learner_folds))
 rows<-x$folds[[1L]]$training_set;H<-prepared$task$vars$history('A',2L)
 A<-prepared$task$vars$A[[2L]];B<-as.matrix(prepared$encoded[rows,c(H,A),drop=FALSE])
 models<-list('Original weights'=readRDS(file.path(base,'ensemble-stability-v7-full-equations',
   folder,'nuisance-fits-fold1.rds'))$beta_fits[[2L]],
   'Landweber weight CV'=readRDS(file.path(cv_root,mechanism,'bridge-fit.rds')))
 if(mechanism=='discrete_dose')models<-c(models,
   list('Fixed 0.01'=readRDS(file.path(weight_root,'weight-0.01','bridge-fit.rds')),
   'Fixed 0.1'=readRDS(file.path(weight_root,'weight-0.1','bridge-fit.rds'))))
 p<-x$population;pt<-study_task(p,g)$task;measured<-p$R3==1
 Vp<-as.matrix(cbind(pt$natural[measured,H,drop=FALSE],p[measured,c('C3_covariate','Y3'),drop=FALSE]))
 truth<-population_truth_functions(p,g)$beta[measured,2L]
 keys<-cell_key(pt$natural[,c(H,A),drop=FALSE]);included<-p$R2[measured]==1
 for(setting in names(models)) {
  m<-models[[setting]]
  stopifnot(length(m$candidate_failures)==0L,identical(m$fold_id,x$learner_folds[[1L]]),
    all(m$weights>=0),abs(sum(m$weights)-1)<1e-10)
  for(name in names(m$candidates)) {
   candidate<-m$candidates[[name]]
   stopifnot(candidate$tuning$link=='inverse_logit')
   pass<-if(name=='landweber')isTRUE(candidate$solver$stopped_by_tolerance) else isTRUE(candidate$solver$converged)
   solvers[[length(solvers)+1L]]<-data.table(mechanism,configuration=setting,candidate=name,
     tolerance_reached=pass,gradient=candidate$solver$coefficient_gradient,
     iterations=if(is.null(candidate$solver$iterations))NA_integer_ else candidate$solver$iterations)
   if(name!='landweber')stopifnot(pass)
   # Landweber also has an explicit iteration regularization limit. Retain
   # tolerance-not-reached cases as flagged exploratory diagnostics; do not
   # claim they passed the stricter final-check tolerance audit.
   if(name%in%c('sieve_md','landweber'))stopifnot(candidate$instrument_spec$type=='cell')
   cv<-candidate$penalty_cv
   if(!is.null(cv))stopifnot(identical(candidate$penalty_ids,prepared$task$id[rows]),
     length(intersect(candidate$penalty_ids,prepared$task$id[x$folds[[1L]]$validation_set]))==0L,
     candidate$tuning[[candidate$penalty_parameter]]==cv$scale[cv$selected])
  }
  gram<-matrix(0,3L,3L)
  for(fold in sort(unique(m$fold_id))) {
   index<-which(m$fold_id==fold);r<-m$cv_residuals[index,,drop=FALSE];n<-nrow(r)
   direct<-matrix(0,3L,3L)
   for(group in base::split(seq_len(n),cell_key(B[index,,drop=FALSE]))) {
    if(length(group)>1L)direct<-direct+
      (tcrossprod(colSums(r[group,,drop=FALSE]))-crossprod(r[group,,drop=FALSE]))/(n*(length(group)-1L))
   }
   stopifnot(max(abs(direct-m$fold_gram[[fold]]))<1e-10)
   gram<-gram+n/m$n_scored*direct
  }
  difference<-max(abs(gram-m$raw_gram))
  stopifnot(difference<1e-10,min(eigen(m$gram,symmetric=TRUE)$values)>-1e-10)
  checked<-0L
  for(scope in seq_len(length(m$fold_penalty_cv)+1L)) {
   cvs<-if(scope==1L)lapply(m$candidates,`[[`,'penalty_cv') else m$fold_penalty_cv[[scope-1L]]
   for(name in names(cvs)) {
    cv<-cvs[[name]];if(is.null(cv))next
    chosen<-which(cv$selected);valid<-is.finite(cv$loss)&!cv$failed
    stopifnot(length(chosen)==1L,chosen>1L,chosen<nrow(cv),all(cv$scale>0),
      valid[chosen],cv$loss[chosen]==min(cv$loss[valid]))
    selections[[length(selections)+1L]]<-data.table(mechanism,configuration=setting,
      training_sample=if(scope==1L)'Final' else paste('Learner',scope-1L),candidate=name,
      parameter=if(name=='landweber')'weight_ridge' else 'lambda',
      selected=as.numeric(cv$scale[chosen]),loss=cv$loss[chosen],
      failed_trials=sum(cv$failed),interior=TRUE)
    checked<-checked+1L
    if(scope==1L)curves[[length(curves)+1L]]<-cbind(as.data.table(cv),
      data.table(mechanism,configuration=setting,candidate=name))
   }
  }
  audits[[length(audits)+1L]]<-data.table(mechanism,configuration=setting,
    raw_gram_max_difference=difference,positive_interior_selections=checked,
    candidate_failures=length(m$candidate_failures),projection_norm=m$psd_projection$adjustment_norm)
  weights[[length(weights)+1L]]<-data.table(mechanism,configuration=setting,
    candidate=names(m$weights),weight=as.numeric(m$weights))
  predictors<-c(list(ensemble=m),m$candidates)
  for(name in names(predictors)) {
   fitted<-as.numeric(predict(predictors[[name]],Vp))
   stopifnot(all(is.finite(fitted)),all(fitted>=1))
   residual<-rep(-1,nrow(p));residual[measured]<-fitted-1
   equation<-data.table(cell=keys,probability=p$probability,residual=residual)
   equation<-equation[,.(probability=sum(probability),error=sum(probability*residual)/sum(probability)),by=cell]
   metrics[[length(metrics)+1L]]<-data.table(mechanism,configuration=setting,candidate=name,
     bridge_rmse=sqrt(sum(p$probability[measured][included]*(fitted[included]-truth[included])^2)/sum(p$probability[measured][included])),
     equation_rmse=sqrt(sum(equation$probability*equation$error^2)),
     equation_maximum_error=max(abs(equation$error)),minimum=min(fitted),maximum=max(fitted),
     application_bound_n=sum(fitted<=-100|fitted>=100))
   points[[length(points)+1L]]<-data.table(mechanism,configuration=setting,candidate=name,
     reference=truth[included],fitted=fitted[included],probability=p$probability[measured][included])
  }
 }
}
for(pair in list(list(metrics,'bridge-comparison.csv'),list(weights,'ensemble-weights.csv'),
 list(audits,'fit-audit.csv'),list(selections,'positive-penalty-audit.csv'),list(curves,'penalty-curves.csv'),list(solvers,'solver-audit.csv')))
 fwrite(rbindlist(pair[[1L]],fill=TRUE),file.path(weight_root,pair[[2L]]))
points<-rbindlist(points);fwrite(points,file.path(weight_root,'bridge-plot-points.csv'))
dir.create(file.path(weight_root,'figures'),showWarnings=FALSE)
for(mechanism_name in c('binary_longitudinal','discrete_dose')) {
 shown<-points[mechanism==mechanism_name];limits<-range(shown$reference,shown$fitted)
 shown[,candidate:=factor(candidate,levels=c('ensemble','sieve_md','landweber','pmmr'),
  labels=c('Ensemble','Sieve minimum distance','Landweber','PMMR'))]
 figure<-ggplot(shown,aes(reference,fitted,size=probability))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey35')+
  geom_point(alpha=.55,colour='#21618b')+facet_grid(configuration~candidate)+
  coord_equal(xlim=limits,ylim=limits)+scale_size_continuous(range=c(.7,3),guide='none')+
  labs(x='Population bridge solution',y='Estimated bridge',
   title=paste(if(mechanism_name=='discrete_dose')'Numerical dose' else 'Binary treatment',
     'n = 4,000, replication 6, first training sample'),
   subtitle='Outcome at time 3; intermediate and final visits observed',
   caption='All predictions retained on equal axes. Fixed-weight checks are exploratory; CV uses only training samples.')+
  theme_bw(base_size=10)
 for(ext in c('png','pdf'))ggsave(file.path(weight_root,'figures',paste0(mechanism_name,'-weight-comparison.',ext)),
   figure,width=14,height=if(mechanism_name=='discrete_dose')14 else 8,dpi=180)
}
cat('Audited six saved models, positive interior selections, shared splits and independent cell U-Grams; no refits.\n')
