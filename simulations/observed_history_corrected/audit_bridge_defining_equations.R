# Population algebra only: check which equations the bridge fitting criterion
# identifies. No datasets, known-function estimators, or zero-penalty fits.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
output<-Sys.getenv('BRIDGE_EQUATION_AUDIT_OUTPUT');stopifnot(nzchar(output))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
summary<-examples<-list()
for(mechanism in c('binary_longitudinal','discrete_dose')) {
  g<-make_mechanism(mechanism,dose_max=3L,visits='missing_both')
  require_valid_mechanism(g)
  p<-enumerate_data(g);prepared<-study_task(p,g);pt<-prepared$task
  reference<-population_truth_functions(p,g);w<-p$probability
  for(s in 1:2) {
    H<-pt$vars$history('A',s);A<-pt$vars$A[[s]]
    B<-as.matrix(pt$natural[,c(H,A),drop=FALSE])
    measured<-p[[paste0('R',s+1L)]]==1
    health<-c(paste0('C',s+1L,'_covariate'),paste0('Y',s+1L))
    V<-as.matrix(cbind(pt$natural[measured,H,drop=FALSE],p[measured,health,drop=FALSE]))
    bkeys<-cell_key(B);vkeys<-cell_key(V)
    blevels<-unique(bkeys);vlevels<-unique(vkeys)
    mass<-as.numeric(rowsum(w,match(bkeys,blevels),reorder=FALSE))
    joint<-as.matrix(Matrix::sparseMatrix(i=match(bkeys[measured],blevels),
      j=match(vkeys,vlevels),x=w[measured],dims=c(length(blevels),length(vlevels))))
    conditional<-joint/mass
    polynomial<-cmbridge::poly_basis(B,degree=1L)
    integrated<-as.matrix(crossprod(polynomial[measured,,drop=FALSE],
      Matrix::sparseMatrix(i=seq_along(vkeys),j=match(vkeys,vlevels),
        x=w[measured],dims=c(length(vkeys),length(vlevels)))))
    right<-as.numeric(crossprod(polynomial,w))
    beta<-reference$beta[which(measured)[match(vlevels,vkeys)],s]
    stopifnot(max(abs(conditional%*%beta-1))<1e-10,
      max(abs(integrated%*%beta-right))<1e-10)
    z<-svd(integrated)
    keep<-z$d>max(z$d)*1e-12
    rowspace<-z$v[,keep,drop=FALSE]
    main_design<-cmbridge::poly_basis(V[match(vlevels,vkeys),,drop=FALSE],degree=1L)
    main_singular<-svd(main_design,nu=0L,nv=0L)$d
    main_jacobian<-integrated%*%((beta-1)*main_design)
    main_jacobian_singular<-svd(main_jacobian,nu=0L,nv=0L)$d
    remaining<-conditional-(conditional%*%rowspace)%*%t(rowspace)
    chosen<-which.max(mass*rowSums(remaining^2))
    direction<-as.numeric(remaining[chosen,])
    stopifnot(max(abs(direction))>1e-8)
    # Stay strictly above one. The alternative remains in the unrestricted
    # inverse-expit joint-category target class on the complete support.
    alternative<-beta+.5*min(beta-1)*direction/max(abs(direction))
    spec<-cmbridge:::.fit_basis_spec(V,'cell_linear',1L,6L)
    coefficients<-numeric(ncol(cmbridge:::.eval_basis_spec(V,spec)))
    coefficients[seq.int(spec$p+2L,length(coefficients))]<-
      -log(alternative[match(spec$levels,vlevels)]-1)
    represented<-1+exp(-as.numeric(cmbridge:::.eval_basis_spec(V,spec)%*%coefficients))
    stopifnot(min(alternative)>1,
      max(abs(represented-alternative[match(vkeys,vlevels)]))<1e-10)
    error<-as.numeric(conditional%*%alternative)-1
    polynomial_error<-max(abs(integrated%*%alternative-right))
    singular<-svd(conditional,nu=0L,nv=0L)$d
    summary[[length(summary)+1L]]<-data.table(mechanism=mechanism,outcome_time=s+1L,
      conditioning_combinations=length(blevels),target_combinations=length(vlevels),
      polynomial_columns=ncol(polynomial),polynomial_operator_rank=sum(keep),
      full_conditional_operator_rank=sum(singular>max(singular)*1e-12),
      landweber_main_target_rank=sum(main_singular>max(main_singular)*1e-12),
      polynomial_moment_jacobian_rank_for_main_target=
        sum(main_jacobian_singular>max(main_jacobian_singular)*1e-12),
      reference_full_equation_max_error=max(abs(conditional%*%beta-1)),
      alternative_polynomial_moment_max_error=polynomial_error,
      alternative_full_equation_rmse=sqrt(sum(mass*error^2)),
      alternative_full_equation_max_error=max(abs(error)),
      alternative_minimum=min(alternative),alternative_maximum=max(alternative),
      target_class_representation_max_error=max(abs(represented-alternative[match(vkeys,vlevels)])))
    examples[[length(examples)+1L]]<-data.table(mechanism=mechanism,outcome_time=s+1L,
      conditioning_cell=blevels,probability=mass,alternative_equation_error=error)
    stopifnot(polynomial_error<1e-9,sqrt(sum(mass*error^2))>1e-4)
  }
}
summary<-rbindlist(summary)
fwrite(summary,file.path(output,'bridge-equation-identification.csv'))
fwrite(rbindlist(examples),file.path(output,'alternative-conditional-equations.csv'))
print(summary)
cat('Different positive functions satisfy all polynomial moments but violate full defining equations.\n')
