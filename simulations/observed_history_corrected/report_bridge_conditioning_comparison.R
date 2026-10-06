# Read-only completion and graphs for the two already saved diagnostic fits.
# This never fits a nuisance model or changes the frozen study.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
library(ggplot2)
root<-Sys.getenv('BRIDGE_CONDITIONING_ROOT')
input<-Sys.getenv('BRIDGE_CONDITIONING_INPUT')
stopifnot(nzchar(root),nzchar(input))
x<-readRDS(input)
g<-make_mechanism(x$design$mechanism,dose_max=3L,visits='missing_both')
p<-x$population;pt<-study_task(p,g)$task
H<-pt$vars$history('A',2L);measured<-p$R3==1
Vp<-as.matrix(cbind(pt$natural[measured,H,drop=FALSE],
  p[measured,c('C3_covariate','Y3'),drop=FALSE]))
baseline<-x$fitted_functions[['sdr/beta1_lambda1_ratio1_m1/1']]$beta[measured,2L]
truth<-population_truth_functions(p,g)$beta[measured,2L]
records<-plots<-weights<-selected_penalties<-solver_records<-list()
for(b in c('polynomial','joint-categories')) {
  path<-file.path(root,b);m<-readRDS(file.path(path,'bridge-fit.rds'))
  stopifnot(length(m$candidate_failures)==0L,
    all(is.finite(m$weights)),all(m$weights>=0),abs(sum(m$weights)-1)<1e-10)
  for(name in names(m$candidates)) {
    solver<-m$candidates[[name]]$solver
    converged<-if(name=='landweber')isTRUE(solver$stopped_by_tolerance) else isTRUE(solver$converged)
    stopifnot(converged,is.finite(solver$coefficient_gradient))
    solver_records[[length(solver_records)+1L]]<-data.table(basis=b,candidate=name,
      recorded_convergence=converged,coefficient_gradient=solver$coefficient_gradient,
      final_loss=solver$final_loss,criterion=if(name=='landweber')'stopped_by_tolerance' else 'converged')
  }
  original_difference<-max(abs(as.numeric(predict(m,Vp))-baseline))
  if(b=='polynomial')stopifnot(original_difference<1e-9)
  checks<-list(m$candidates$sieve_md$penalty_cv)
  checks<-c(checks,unlist(m$fold_penalty_cv,recursive=FALSE))
  checked<-0L
  for(cv in checks) {
    if(is.null(cv)||!nrow(cv)||!all(c('scale','loss','selected')%in%names(cv)))next
    chosen<-cv[cv$selected,,drop=FALSE]
    valid<-cv[is.finite(cv$loss)&!cv$failed,,drop=FALSE]
    stopifnot(nrow(chosen)==1L,is.finite(chosen$scale),chosen$scale>0,
      !chosen$failed,chosen$loss==min(valid$loss),
      chosen$scale>min(valid$scale),chosen$scale<max(valid$scale))
    checked<-checked+1L
    selected_penalties[[length(selected_penalties)+1L]]<-data.table(basis=b,
      fit=if(checked==1L)'final candidate' else paste('learner training sample',checked-1L),
      scale=chosen$scale,loss=chosen$loss,failed_trial_records=sum(cv$failed),
      grid_points=nrow(cv),interior=TRUE)
  }
  stopifnot(checked==4L)
  points<-fread(file.path(path,'bridge-plot-points.csv'))
  for(name in unique(points$candidate)) {
    fitted<-as.numeric(predict(if(name=='ensemble')m else m$candidates[[name]],Vp))
    saved<-points[candidate==name]
    selected<-p$R2[measured]==1
    stopifnot(max(abs(saved$fitted-fitted[selected]))<1e-10,
      max(abs(saved$reference-truth[selected]))<1e-10)
  }
  records[[b]]<-data.table(basis=b,baseline_prediction_max_difference=original_difference,
    candidate_failure_records=length(m$candidate_failures),penalty_cv_records_checked=checked,
    recovered_reporting_tail=TRUE)
  points[,basis_label:=if(b=='polynomial')'Original conditioning basis' else 'Joint-category conditioning basis']
  plots[[b]]<-points
  w<-fread(file.path(path,'ensemble-weights.csv'));w[,basis:=b];weights[[b]]<-w
}
fwrite(rbindlist(records),file.path(root,'saved-fit-audit.csv'))
fwrite(rbindlist(weights),file.path(root,'ensemble-weights.csv'))
fwrite(rbindlist(selected_penalties),file.path(root,'selected-penalty-audit.csv'))
fwrite(rbindlist(solver_records),file.path(root,'solver-audit.csv'))
metrics<-rbindlist(lapply(c('polynomial','joint-categories'),function(b)
  fread(file.path(root,b,'population-bridge-checks.csv'))))
fwrite(metrics,file.path(root,'population-bridge-checks.csv'))
points<-rbindlist(plots)
points[,candidate:=factor(candidate,levels=c('ensemble','sieve_md','landweber','pmmr'),
  labels=c('Ensemble','Sieve minimum distance','Landweber','PMMR'))]
points[,basis_label:=factor(basis_label,
  levels=c('Original conditioning basis','Joint-category conditioning basis'))]
limits<-range(points$reference,points$fitted)
figure<-ggplot(points,aes(reference,fitted,size=probability))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+
  geom_point(alpha=.55,colour='#215d80')+
  facet_grid(basis_label~candidate)+coord_equal(xlim=limits,ylim=limits)+
  scale_size_continuous(range=c(.6,3),guide='none')+
  labs(x='Population bridge solution',y='Estimated bridge',
    title='Binary treatment, n = 4,000, replication 2, first training sample',
    subtitle='Outcome at time 3; R2 = 1 and R3 = 1',
    caption='All combinations are shown on the same equal axes. Point size reflects population probability. No refitting in this report.')+
  theme_bw(base_size=11)+theme(strip.text.y=element_text(size=9))
dir.create(file.path(root,'figures'),showWarnings=FALSE)
for(ext in c('png','pdf'))ggsave(file.path(root,'figures',paste0('all-bridge-candidates.',ext)),
  figure,width=14,height=8.5,dpi=180)
writeLines(c('The initial fits completed and were saved successfully.',
  'Both original script invocations then failed at the final reporting expression because of a missing parenthesis.',
  'The two attempt logs are preserved. After repairing that expression, the retry guard refused to overwrite the saved models.',
  'This separate read-only script audited the existing models and completed the comparison; it did not refit them.',
  'The polynomial fit reproduces the original saved ensemble predictions within 1e-9.',
  'No final candidate failures, application-bound hits, or invalid selected penalties.'),file.path(root,'AUDIT.txt'))
cat('Audited two saved models and completed graphs without refitting.\n')
