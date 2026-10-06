# Partition saved time-3 population equations by observed-history status.
# All estimated functions remain fixed. No nuisance parameters are refitted.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
input<-Sys.getenv('DIAGNOSTIC_INPUT');output<-Sys.getenv('DIAGNOSTIC_OUTPUT')
stopifnot(nzchar(input),nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
records<-list()
for(path in list.files(file.path(input,'jobs'),'\\.rds$',full.names=TRUE)) {
  x<-readRDS(path);p<-x$population
  if(is.null(p)||!length(x$fitted_functions))next
  g<-make_mechanism(x$design$mechanism,dose_max=3L,visits='missing_both')
  pt<-study_task(p,g)$task;truth<-population_truth_functions(p,g)
  H<-pt$vars$history('A',2L);A<-pt$vars$A[[2L]]
  B<-cell_key(pt$natural[,c(H,A),drop=FALSE])
  measured<-p$R3==1;w<-p$probability
  V<-cell_key(cbind(pt$natural[measured,H,drop=FALSE],
    p[measured,c('C3_covariate','Y3'),drop=FALSE]))
  omega<-apply(truth$ratios,1,prod)
  for(j in seq_along(x$folds)) {
    nu<-x$fitted_functions[[paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')]]
    if(is.null(nu))next
    remainder<-numeric(nrow(p))
    remainder[measured]<--(nu$lambda[measured,2L]-truth$lambda[measured,2L])*
      (nu$beta[measured,2L]-truth$beta[measured,2L])
    residual<-rep(-1,nrow(p));residual[measured]<-nu$beta[measured,2L]-1
    adjoint<-data.table(R2=p$R2[measured],cell=V,probability=w[measured],
      residual=nu$lambda[measured,2L]-omega[measured]*p$Y3[measured])
    adjoint<-adjoint[,.(probability=sum(probability),
      equation_error=sum(probability*residual)/sum(probability)),by=.(R2,cell)]
    adjoint<-adjoint[,.(adjoint_equation_rmse=
      sqrt(sum(probability*equation_error^2)/sum(probability))),by=R2]
    bridge<-data.table(R2=p$R2,cell=B,probability=w,residual=residual)
    bridge<-bridge[,.(probability=sum(probability),
      equation_error=sum(probability*residual)/sum(probability)),by=.(R2,cell)]
    bridge<-bridge[,.(bridge_equation_rmse=
      sqrt(sum(probability*equation_error^2)/sum(probability))),by=R2]
    contribution<-data.table(R2=p$R2,probability=w,measured=measured,remainder=remainder,
      squared_beta_error=(nu$beta[,2L]-truth$beta[,2L])^2,
      squared_lambda_error=(nu$lambda[,2L]-truth$lambda[,2L])^2)
    contribution<-contribution[,.(population_probability=sum(probability),
      measured_probability=sum(probability*measured),
      bridge_remainder_contribution=sum(probability*remainder),
      beta_rmse_measured=sqrt(sum(probability[measured]*squared_beta_error[measured])/
        sum(probability[measured])),
      lambda_rmse_to_selected_solution_measured=
        sqrt(sum(probability[measured]*squared_lambda_error[measured])/sum(probability[measured]))),by=R2]
    contribution<-merge(merge(contribution,bridge,by='R2'),adjoint,by='R2')
    reference<-x$diagnostics[estimator=='sdr' & fold==j & horizon==2]
    stopifnot(nrow(reference)==1L,
      abs(sum(contribution$bridge_remainder_contribution)-reference$bridge_remainder)<1e-8,
      abs(sum(contribution$population_probability)-1)<1e-12)
    contribution[,`:=`(mechanism=x$design$mechanism,n=x$design$n,
      replicate=x$design$replicate,seed=x$design$seed,fold=j,
      fold_weight=length(x$folds[[j]]$validation_set)/x$design$n)]
    records[[length(records)+1L]]<-contribution
  }
}
records<-rbindlist(records)
fwrite(records,file.path(output,'time3-R2-contributions-by-fold.csv'))
columns<-setdiff(names(records),c('mechanism','n','replicate','seed','fold','fold_weight','R2'))
combined<-records[,lapply(.SD,function(z)sum(z*fold_weight)),
  by=.(mechanism,n,replicate,seed,R2),.SDcols=columns]
fwrite(combined,file.path(output,'time3-R2-contributions-by-dataset.csv'))
cat('Time-3 subgroup contributions add exactly to each saved bridge remainder.\n')
