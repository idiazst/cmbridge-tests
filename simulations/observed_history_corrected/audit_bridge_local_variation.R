# Algebraic local variation implied by the full conditional moments.
# No oracle fit or known-function performance simulation is run.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
base<-Sys.getenv('BRIDGE_VARIATION_BASE');output<-Sys.getenv('BRIDGE_VARIATION_OUTPUT')
stopifnot(nzchar(base),nzchar(output));records<-list()
for(mechanism in c('binary_longitudinal','discrete_dose')) {
 seed<-if(mechanism=='binary_longitudinal')5103006L else 5203006L
 folder<-paste0(mechanism,'-n4000-seed',seed)
 x<-readRDS(file.path(base,'ensemble-stability-v6',folder,'result.rds'))
 model<-readRDS(file.path(base,'ensemble-stability-v8-weight-cv',folder,'nuisance-fits-fold2.rds'))$beta_fits[[2L]]$candidates$landweber
 p<-x$population;g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
 pt<-study_task(p,g)$task;H<-pt$vars$history('A',2L);A<-pt$vars$A[[2L]]
 B<-as.matrix(pt$natural[,c(H,A),drop=FALSE]);observed<-p$R3==1
 V<-as.matrix(cbind(pt$natural[observed,H,drop=FALSE],p[observed,c('C3_covariate','Y3'),drop=FALSE]))
 beta<-population_truth_functions(p,g)$beta[observed,2L]
 features<-cmbridge:::.eval_basis_spec(V,model$target_spec)
 # Remove exact algebraic coefficient redundancies ONLY in this matrix audit.
 # This does not modify fitted models, predictors, basis classes or penalties.
 decomposition<-svd(features*sqrt(p$probability[observed]))
 rank<-sum(decomposition$d>max(decomposition$d)*1e-10)
 reduced<-features%*%decomposition$v[,seq_len(rank),drop=FALSE]
 derivative<-matrix(0,nrow(p),rank);derivative[observed,]<--(beta-1)*reduced
 residual<-rep(-1,nrow(p));residual[observed]<-beta-1
 index<-match(cell_key(B),unique(cell_key(B)))
 probability<-as.numeric(rowsum(p$probability,index,reorder=FALSE))
 sigma<-as.numeric(rowsum(p$probability*residual^2,index,reorder=FALSE))/probability
 jacobian<-rowsum(derivative*p$probability,index,reorder=FALSE)/probability
 mean_residual<-as.numeric(rowsum(p$probability*residual,index,reorder=FALSE))/probability
 stopifnot(max(abs(mean_residual))<1e-10,all(sigma>0))
 information<-crossprod(jacobian*sqrt(probability/sigma))
 eigen<-eigen(information,symmetric=TRUE)
 stopifnot(min(eigen$values)>max(eigen$values)*1e-12)
 sample_n<-length(x$folds[[2L]]$training_set)
 covariance<-eigen$vectors%*%diag(1/eigen$values,rank)%*%t(eigen$vectors)/sample_n
 variance<-rowSums((derivative%*%covariance)*derivative)
 included<-p$R3==1 & p$R2==1
 expected_prediction_rmse<-sqrt(sum(p$probability[included]*variance[included])/sum(p$probability[included]))
 records[[length(records)+1L]]<-data.table(mechanism,n=4000L,training_rows=sample_n,
   independent_target_parameters=rank,conditioning_combinations=length(probability),
   bridge_equation_maximum_error=max(abs(mean_residual)),
   local_information_minimum_eigenvalue=min(eigen$values),
   local_information_condition_number=max(eigen$values)/min(eigen$values),
   local_optimally_weighted_prediction_rmse=expected_prediction_rmse)
}
fwrite(rbindlist(records),file.path(output,'local-moment-variation.csv'));print(rbindlist(records))
