source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
output <- Sys.getenv('SIM_OUTPUT'); rows <- list()
for (name in c('binary_longitudinal','discrete_dose')) {
 g<-make_mechanism(name,dose_max=3L,visits='missing_both');p<-enumerate_data(g);pt<-study_task(p,g)$task;tru<-population_truth_functions(p,g)
 for(s in 1:2) {
  H<-pt$vars$history('A',s);obs<-p[[paste0('R',s+1L)]]==1
  x<-as.matrix(cbind(pt$natural[obs,H,drop=FALSE],p[obs,study_arguments(g)$health[[s]],drop=FALSE]))
  # Inverse measurement probability is a valid solution at every history,
  # including histories with deterministic treatment and nonunique bridges.
  current_y<-p[[paste0('Y',s)]][obs];current_y[is.na(current_y)]<-0
  eta<-.25-.65*p[[paste0('Y',s+1L)]][obs]-.45*p[[paste0('C',s+1L,'_covariate')]][obs]+.15*current_y+.1*p$C1_baseline[obs]
  beta<-1+exp(-eta);R<-p[[paste0('R',s)]][obs]
  stopifnot(max(abs(beta[R==1]-tru$beta[obs,s][R==1]))<1e-12)
  # Construct the basis from just one observed training combination; the
  # unseen combinations must still have the correct fixed main-effect span.
  for(a in c('sieve_md','landweber')) {
   control<-cm_library(bridge=TRUE)[[a]]$control
   spec<-cmbridge:::.fit_basis_spec(x[1,,drop=FALSE],control$target_basis,if (is.null(control$target_degree)) 3L else control$target_degree)
   design<-cmbridge:::.eval_basis_spec(x,spec);dec<-svd(design);keep<-dec$d>max(dec$d)*1e-10
   coefficients<-as.numeric(dec$v[,keep,drop=FALSE]%*%(as.numeric(crossprod(dec$u[,keep,drop=FALSE],eta))/dec$d[keep]))
   represented<-1+exp(-as.numeric(design%*%coefficients));error<-max(abs(represented-beta));stopifnot(error<1e-10)
   for(state in unique(R)) rows[[length(rows)+1L]]<-data.table(mechanism=name,horizon=s,candidate=a,current_R=state,reachable_combinations=uniqueN(cell_key(x[R==state,,drop=FALSE])),basis=spec$type,maximum_representation_error=max(abs(represented[R==state]-beta[R==state])),basis_built_from_one_combination=TRUE)
  }
 }
}
fwrite(rbindlist(rows),file.path(output,'bridge-function-class.csv'));print(rbindlist(rows))
cat('Both classes contain a valid bridge solution on the full reachable support. No truth used in fitting.\n')
