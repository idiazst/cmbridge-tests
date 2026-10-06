# Algebraic expectation of the existing training objective at a valid bridge.
# No nuisance model or known-function performance simulation is fitted.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
reference<-Sys.getenv('BRIDGE_GRADIENT_REFERENCE');output<-Sys.getenv('BRIDGE_GRADIENT_OUTPUT')
stopifnot(nzchar(reference),nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
records<-list()
for(mechanism in c('binary_longitudinal','discrete_dose')) {
  seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
  x<-readRDS(file.path(reference,paste0(mechanism,'-n4000-seed',seed),'result.rds'))
  g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both');require_valid_mechanism(g)
  d<-x$data;p<-x$population
  task<-study_task(d,g,folds=x$folds,learner_folds=x$learner_folds,learner_groups=3L)$task
  pt<-study_task(p,g)$task
  H<-pt$vars$history('A',2L);A<-pt$vars$A[[2L]]
  Bpop<-as.matrix(pt$natural[,c(H,A),drop=FALSE])
  measured<-p$R3==1;beta<-numeric(nrow(p))
  current_y<-p$Y2[measured];current_y[is.na(current_y)]<-0
  eta<-.25-.65*p$Y3[measured]-.45*p$C3_covariate[measured]+.15*current_y+.1*p$C1_baseline[measured]
  beta[measured]<-1+exp(-eta)
  residual<-rep(-1,nrow(p));residual[measured]<-beta[measured]-1
  # At a valid bridge, E[residual | conditioning values] = 0.
  # Increasing the unpenalized intercept reduces beta by beta-1.
  # Conditional covariance(residual, d residual/d intercept) is
  # -E[R3*(beta-1)^2 | conditioning values].
  moments<-data.table(cell=cell_key(Bpop),probability=p$probability,
    residual=residual,squared=residual^2,gap_squared=as.numeric(measured)*(beta-1)^2)
  moments<-moments[,.(probability=sum(probability),
    conditional_residual=sum(probability*residual)/sum(probability),
    residual_variance=sum(probability*squared)/sum(probability),
    measured_squared_gap=sum(probability*gap_squared)/sum(probability)),by=cell]
  equation_error<-max(abs(moments$conditional_residual));stopifnot(equation_error<1e-10)
  training<-x$folds[[1L]]$training_set
  B<-as.matrix(task$natural[training,c(H,A),drop=FALSE]);n<-nrow(B)
  M<-d$R3[training];active<-M==1
  actual_beta<-numeric(n);current_y<-d$Y2[training][active];current_y[is.na(current_y)]<-0
  actual_eta<-.25-.65*d$Y3[training][active]-.45*d$C3_covariate[training][active]+
    .15*current_y+.1*d$C1_baseline[training][active]
  actual_beta[active]<-1+exp(-actual_eta)
  actual_residual<-M*actual_beta-1
  actual_derivative<--M*(actual_beta-1)
  matched<-match(cell_key(B),moments$cell);stopifnot(!anyNA(matched))
  settings<-data.table(basis=c('poly','cell','cell','cell'),
    ridge=c(1e-8,1e-8,.01,.1))
  for(i in seq_len(nrow(settings))) {
    setting<-settings[i]
    spec<-cmbridge:::.fit_basis_spec(B,setting$basis,1L,6L)
    critic<-cmbridge:::.moment_basis_spec(B,spec,setting$ridge)
    influence<-rowSums((critic$Q%*%critic$W)*critic$Q)
    stopifnot(all(is.finite(influence)),min(influence)>-1e-8)
    floor<-sum(influence*moments$residual_variance[matched])/n^2
    gradient<--2*sum(influence*moments$measured_squared_gap[matched])/n^2
    empirical_moment<-as.numeric(crossprod(critic$Q,actual_residual))/n
    empirical_derivative<-as.numeric(crossprod(critic$Q,actual_derivative))/n
    empirical_gradient<-2*as.numeric(crossprod(empirical_moment,critic$W%*%empirical_derivative))
    empirical_loss<-as.numeric(crossprod(empirical_moment,critic$W%*%empirical_moment))
    stopifnot(floor>0,gradient<0)
    records[[length(records)+1L]]<-data.table(mechanism=mechanism,n=4000L,seed=seed,
      outer_training_sample=1L,outcome_time=3L,training_rows=n,
      conditioning_basis=setting$basis,conditioning_weight_ridge=setting$ridge,
      conditioning_combinations=uniqueN(cell_key(B)),
      valid_bridge_maximum_conditional_error=equation_error,
      expected_training_moment_loss_at_valid_bridge=floor,
      expected_intercept_derivative_at_valid_bridge=gradient,
      actual_training_moment_loss_at_valid_bridge=empirical_loss,
      actual_intercept_derivative_at_valid_bridge=empirical_gradient,
      maximum_training_row_weight=max(influence))
  }
}
fwrite(rbindlist(records),file.path(output,'expected-training-gradient.csv'))
writeLines(c('Analytic expectation conditional on the existing training conditioning values.',
  'No nuisance fit, new simulated dataset, known-function performance study, or zero-penalty simulation.',
  'The inverse-expit bridge used in the calculation is a valid full-support solution; its conditional equation is independently checked.',
  'At that solution, the expected squared empirical moment objective contains a function-dependent variance term.',
  'The negative intercept derivative favors decreasing the bridge, although all conditional equations hold at the reference solution.',
  'The link ridge does not penalize the intercept, so its derivative in this direction is zero for every positive chosen lambda.',
  'Conditioning-weight stabilization reduces this finite-sample contribution without deleting equations.',
  'This does not prove that the actual optimizer or coverage errors are fully explained by this one direction.'),
  file.path(output,'NOTES.txt'))
print(rbindlist(records))
