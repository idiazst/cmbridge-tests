# Summarize one completed local runtime diagnosis, excluded from study counts.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
root<-Sys.getenv('RUNTIME_PROFILE_OUTPUT');stopifnot(nzchar(root))
x<-readRDS(file.path(root,'result.rds'))
stopifnot(nrow(x$results)==4L,nrow(x$errors)==0L,x$design$n==4000L,
  x$design$seed==5103002L,x$design$mechanism=='binary_longitudinal')
selected<-x$tuning[selected==TRUE]
stopifnot(nrow(selected)>0L,all(selected$scale>0),all(is.finite(selected$loss)),
  all(!selected$failed),all(startsWith(selected$penalty_status,'interior')))
stopifnot(max(abs(x$results$mean_total_eif-x$results$mean_sequential_eif-
  x$results$mean_bridge_adjoint_eif))<1e-10,
  max(abs(x$diagnostics$bridge_identity_error))<1e-8)
for(name in c('results','diagnostics','tuning','ensemble_weights',
             'candidate_failures','penalty_trial_failures','errors')) {
  value<-x[[name]]
  if(!is.null(value)&&ncol(value))fwrite(value,file.path(root,paste0(name,'.csv')))
}
fwrite(selected,file.path(root,'selected-penalties.csv'))
population<-x$diagnostics[,lapply(.SD,function(z)sum(z*fold_weight)),
  by=.(estimator,horizon),.SDcols=c('population_bias','sequential_remainder','bridge_remainder')]
fwrite(population,file.path(root,'population-by-outcome-time.csv'))
writeLines(unique(x$job_warnings),file.path(root,'warnings.txt'))
g<-make_mechanism(x$design$mechanism,dose_max=3L,visits='missing_both')
p<-x$population;truth<-population_truth_functions(p,g)
points<-metrics<-list()
for(j in seq_along(x$folds))for(s in 1:2) {
  nu<-x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')]]
  observed<-p[[paste0('R',s+1L)]]==1
  if(s==2L)observed<-observed & p$R2==1
  points[[length(points)+1L]]<-data.table(fold=factor(j),outcome_time=paste('Outcome at time',s+1L),
    reference=truth$beta[observed,s],fitted=nu$beta[observed,s])
  metrics[[length(metrics)+1L]]<-data.table(fold=j,outcome_time=s+1L,
    rmse=sqrt(sum(p$probability[observed]*(nu$beta[observed,s]-truth$beta[observed,s])^2)/
      sum(p$probability[observed])),
    application_bound_fraction=mean(nu$beta[observed,s]<=-100 | nu$beta[observed,s]>=100))
}
points<-rbindlist(points);metrics<-rbindlist(metrics)
fwrite(metrics,file.path(root,'bridge-rmse.csv'))
limits<-range(points$reference,points$fitted)
plot<-ggplot(points,aes(reference,fitted))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.3,size=.75,colour='#2876B2')+
  facet_grid(outcome_time~fold)+coord_equal(xlim=limits,ylim=limits)+
  labs(x='Valid bridge solution',y='Fitted ensemble bridge',
    title='Independent local runtime diagnosis: binary treatment, n = 4,000',
    subtitle='Columns show the three training samples',
    caption='All measured predictions shown; time 3 additionally requires R2 = 1. Equal axes; no cropped outliers.')+
  theme_bw(base_size=11)
dir.create(file.path(root,'figures'),showWarnings=FALSE)
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('bridge-ensemble.',ext)),
  plot,width=10,height=7,dpi=180)
equations<-fread(file.path(root,'equations','adjoint-conditional-equations.csv'))
equations[,`:=`(fold=factor(fold),outcome_time=paste('Outcome at time',horizon+1L),
  measured_training_group=factor(ifelse(measured_training_n==0,'0',
    ifelse(measured_training_n<=2,'1 or 2','3 or more')),levels=c('0','1 or 2','3 or more')))]
limits<-range(equations$conditional_true_loading,equations$conditional_fitted_adjoint)
plot<-ggplot(equations,aes(conditional_true_loading,conditional_fitted_adjoint,
  colour=measured_training_group))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.65,size=1.1)+facet_grid(outcome_time~fold)+
  coord_equal(xlim=limits,ylim=limits)+
  labs(x='Right side of the population adjoint equation',
    y='Conditional mean of the fitted adjoint',colour='Measured training rows',
    title='Adjoint defining equations from the saved runtime diagnosis',
    subtitle='Columns show the three training samples',
    caption='Every positive-probability conditioning combination is shown. Equal axes; no cropped outliers.')+
  theme_bw(base_size=11)+theme(legend.position='bottom')
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('adjoint-equations.',ext)),
  plot,width=10,height=7,dpi=180)
cat('Runtime diagnosis audited; all selected penalties interior; saved defining identity verified.\n')
print(metrics);print(population)
