# Diagnostic calculations on already fitted functions. No nuisance fitting.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
input<-Sys.getenv('DIAGNOSTIC_INPUT')
output<-Sys.getenv('DIAGNOSTIC_OUTPUT')
stopifnot(nzchar(input),nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
summaries<-adjoint_cells<-bridge_cells<-list()
projection_error<-function(operator,loading,probability) {
  scale<-sqrt(probability/sum(probability))
  z<-svd(operator*scale)
  keep<-z$d>max(z$d)*1e-12
  coefficients<-as.numeric(z$v[,keep,drop=FALSE]%*%
    (crossprod(z$u[,keep,drop=FALSE],loading*scale)/z$d[keep]))
  residual<-as.numeric(operator%*%coefficients)-loading
  c(rmse=sqrt(sum(scale^2*residual^2)),maximum=max(abs(residual)),rank=sum(keep),
    condition_number=max(z$d)/min(z$d[keep]))
}
paths<-list.files(file.path(input,'jobs'),'\\.rds$',full.names=TRUE)
stopifnot(length(paths)>0L)
for(path in paths) {
  x<-readRDS(path)
  if(is.null(x$fitted_functions)||!length(x$fitted_functions))next
  config<-x$design
  set.seed(config$seed)
  g<-make_mechanism(config$mechanism,dose_max=3L,visits='missing_both')
  require_valid_mechanism(g)
  stopifnot(max(abs(exact_truth(g)-backward_truth(g)))<1e-12)
  d<-draw_data(config$n,g)
  prepared<-study_task(d,g,folds=3L,learner_groups=3L,balance_measurement=TRUE)
  stopifnot(identical(prepared$task$folds,x$folds),
    identical(prepared$task$learner_folds,x$learner_folds))
  p<-x$population;pt<-study_task(p,g)$task;tru<-population_truth_functions(p,g)
  task<-prepared$task;w<-p$probability
  for(j in seq_along(x$folds)) {
    key<-paste('sdr','beta1_lambda1_ratio1_m1',j,sep='/')
    nu<-x$fitted_functions[[key]]
    if(is.null(nu))next
    train<-x$folds[[j]]$training_set
    weight<-length(x$folds[[j]]$validation_set)/config$n
    for(s in 1:2) {
      H<-pt$vars$history('A',s);A<-pt$vars$A[[s]]
      observed<-p[[paste0('R',s+1L)]]==1
      train_measured<-train[d[[paste0('R',s+1L)]][train]==1]
      B<-cell_key(pt$natural[,c(H,A),drop=FALSE])
      Bd<-cell_key(task$natural[train,c(H,A),drop=FALSE])
      Bd_measured<-cell_key(task$natural[train_measured,c(H,A),drop=FALSE])
      health<-c(paste0('C',s+1L,'_covariate'),paste0('Y',s+1L))
      V<-cell_key(cbind(pt$natural[observed,H,drop=FALSE],p[observed,health,drop=FALSE]))
      Vd<-cell_key(cbind(task$natural[train_measured,H,drop=FALSE],d[train_measured,health,drop=FALSE]))
      omega_true<-apply(tru$ratios[,seq_len(s),drop=FALSE],1,prod)
      omega_fit<-apply(nu$ratios[,seq_len(s),drop=FALSE],1,prod)
      Y<-p[[paste0('Y',s+1L)]][observed]
      a<-data.table(cell=V,probability=w[observed],
        fitted_adjoint=nu$lambda[observed,s],true_loading=omega_true[observed]*Y,
        fitted_loading=omega_fit[observed]*Y,
        error_true=nu$lambda[observed,s]-omega_true[observed]*Y,
        error_fitted=nu$lambda[observed,s]-omega_fit[observed]*Y,
        ratio_difference=(omega_fit[observed]-omega_true[observed])*Y)
      a<-a[,.(probability=sum(probability),
        conditional_fitted_adjoint=sum(probability*fitted_adjoint)/sum(probability),
        conditional_true_loading=sum(probability*true_loading)/sum(probability),
        conditional_fitted_loading=sum(probability*fitted_loading)/sum(probability),
        equation_error=sum(probability*error_true)/sum(probability),
        error_using_fitted_ratios=sum(probability*error_fitted)/sum(probability),
        difference_from_ratio_estimation=sum(probability*ratio_difference)/sum(probability)),by=cell]
      counts<-table(Vd);a[,measured_training_n:=as.integer(counts[match(cell,names(counts))])]
      a[is.na(measured_training_n),measured_training_n:=0L]
      stopifnot(max(abs(a$equation_error-a$error_using_fitted_ratios-a$difference_from_ratio_estimation))<1e-10)
      measured_mass<-sum(a$probability)
      residual<-rep(-1,nrow(p));residual[observed]<-nu$beta[observed,s]-1
      product<-numeric(nrow(p));product[observed]<--w[observed]*
        (nu$lambda[observed,s]-tru$lambda[observed,s])*(nu$beta[observed,s]-tru$beta[observed,s])
      b<-data.table(cell=B,probability=w,measured_probability=w*observed,
        residual=residual,remainder=product,
        lambda=nu$lambda[,s],lambda_reference=tru$lambda[,s])
      b<-b[,.(probability=sum(probability),measured_probability=sum(measured_probability),
        equation_error=sum(probability*residual)/sum(probability),
        bridge_remainder=sum(remainder),adjoint=sum(probability*lambda)/sum(probability),
        selected_true_adjoint=sum(probability*lambda_reference)/sum(probability)),by=cell]
      counts<-table(Bd);b[,training_n:=as.integer(counts[match(cell,names(counts))])]
      b[is.na(training_n),training_n:=0L]
      counts<-table(Bd_measured)
      b[,measured_training_n:=as.integer(counts[match(cell,names(counts))])]
      b[is.na(measured_training_n),measured_training_n:=0L]
      # Algebraic representation audit, not a zero-penalty data fit. Allow all
      # real coefficients and test whether any adjoint solution satisfies the
      # population equations, including the implemented unseen-value rule.
      operator<-as.matrix(Matrix::sparseMatrix(i=match(V,a$cell),j=match(B[observed],b$cell),
        x=w[observed],dims=c(nrow(a),nrow(b))))/a$probability
      full_class<-projection_error(operator,a$conditional_true_loading,a$probability)
      observed_levels<-sort(unique(Bd_measured))
      index<-match(b$cell,observed_levels)
      continuation<-matrix(0,nrow(b),length(observed_levels))
      seen<-which(!is.na(index))
      continuation[cbind(seen,index[seen])]<-1
      if(anyNA(index))continuation[is.na(index),]<-1/length(observed_levels)
      realized_class<-projection_error(operator%*%continuation,
        a$conditional_true_loading,a$probability)
      free_intercept<-matrix(0,nrow(b),length(observed_levels))
      free_intercept[cbind(seen,index[seen])]<-1
      free_intercept_class<-projection_error(operator%*%cbind(1,free_intercept),
        a$conditional_true_loading,a$probability)
      dg<-x$diagnostics[estimator=='sdr' & fold==j & horizon==s]
      stopifnot(nrow(dg)==1L,abs(sum(b$bridge_remainder)-dg$bridge_remainder)<1e-8)
      summary<-data.table(mechanism=config$mechanism,n=config$n,replicate=config$replicate,
        seed=config$seed,fold=j,horizon=s,fold_weight=weight,
        bridge_remainder=sum(b$bridge_remainder),
        bridge_equation_rmse=sqrt(sum(b$probability*b$equation_error^2)),
        adjoint_equation_rmse=sqrt(sum(a$probability*a$equation_error^2)/measured_mass),
        adjoint_error_using_fitted_ratios_rmse=sqrt(sum(a$probability*a$error_using_fitted_ratios^2)/measured_mass),
        adjoint_difference_from_ratio_estimation_rmse=sqrt(sum(a$probability*a$difference_from_ratio_estimation^2)/measured_mass),
        measured_probability_in_unseen_combinations=sum(a[measured_training_n==0L]$probability)/measured_mass,
        measured_probability_in_combinations_with_at_most_two_rows=sum(a[measured_training_n<=2L]$probability)/measured_mass,
        probability_in_unseen_treatment_history_combinations=sum(b[training_n==0L]$probability),
        remainder_from_unseen_treatment_history_combinations=sum(b[training_n==0L]$bridge_remainder),
        remainder_from_combinations_with_at_most_two_rows=sum(b[training_n<=2L]$bridge_remainder),
        measured_probability_in_unseen_adjoint_target_combinations=
          sum(b[measured_training_n==0L]$measured_probability)/measured_mass,
        remainder_from_unseen_adjoint_target_combinations=sum(b[measured_training_n==0L]$bridge_remainder),
        remainder_from_adjoint_target_combinations_with_at_most_two_rows=
          sum(b[measured_training_n<=2L]$bridge_remainder),
        full_support_adjoint_class_minimum_equation_rmse=unname(full_class['rmse']),
        full_support_adjoint_class_maximum_equation_error=unname(full_class['maximum']),
        observed_dictionary_adjoint_class_minimum_equation_rmse=unname(realized_class['rmse']),
        observed_dictionary_adjoint_class_maximum_equation_error=unname(realized_class['maximum']),
        observed_dictionary_with_free_intercept_minimum_equation_rmse=unname(free_intercept_class['rmse']),
        full_support_adjoint_operator_condition_number=unname(full_class['condition_number']),
        observed_dictionary_adjoint_operator_condition_number=unname(realized_class['condition_number']),
        maximum_abs_fitted_adjoint=max(abs(nu$lambda[,s])))
      summaries[[length(summaries)+1L]]<-summary
      for(table in list(a,b))table[,`:=`(mechanism=config$mechanism,n=config$n,replicate=config$replicate,
        seed=config$seed,fold=j,horizon=s,fold_weight=weight)]
      adjoint_cells[[length(adjoint_cells)+1L]]<-a
      bridge_cells[[length(bridge_cells)+1L]]<-b
    }
  }
}
summary<-rbindlist(summaries)
stopifnot(nrow(summary)>0L)
fwrite(summary,file.path(output,'equation-diagnostics-by-fold.csv'))
variables<-setdiff(names(summary),c('mechanism','n','replicate','seed','fold','horizon','fold_weight'))
combined<-summary[,lapply(.SD,function(z)sum(z*fold_weight)),
  by=.(mechanism,n,replicate,seed,horizon),.SDcols=variables]
fwrite(combined,file.path(output,'equation-diagnostics-by-dataset.csv'))
fwrite(rbindlist(adjoint_cells),file.path(output,'adjoint-conditional-equations.csv'))
fwrite(rbindlist(bridge_cells),file.path(output,'bridge-conditional-equations.csv'))
writeLines(c('Diagnostic evaluation of saved fitted functions; no nuisance fits were run.',
  'Fitted-ratio equations use population predictions of the final outer ratio fits.',
  'Actual adjoint training labels used nested out-of-sample ratio predictions; these are not asserted equal.',
  'Remainders use a selected valid true solution only to verify the defining identity.',
  'Adjoint target counts use measured training rows because fit_adjoint fits only on those rows.',
  'Conditioning-value counts and treatment/history target counts refer to different sets of variables.',
  'Class audits solve population linear equations algebraically, never fit nuisance functions to the sample.',
  'The observed dictionary audit includes the implemented mean continuation for unseen target values.',
  'It allows any valid adjoint solution, rather than requiring equality to a chosen reference solution.',
  'The maximum adjoint size is averaged across folds in the dataset table; see individual fold rows for each maximum.'),
  file.path(output,'NOTES.txt'))
print(combined[,.(mechanism,n,replicate,horizon,bridge_remainder,
  adjoint_equation_rmse,adjoint_error_using_fitted_ratios_rmse,
  adjoint_difference_from_ratio_estimation_rmse,
  measured_probability_in_unseen_combinations)])
