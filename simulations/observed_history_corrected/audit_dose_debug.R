# Exact full-history audit for the DGP used in the component check.
source(file.path(Sys.getenv('STUDY_SOURCE'),'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
g <- make_mechanism('discrete_dose',visits='scheduled',dose_max=3L)
require_valid_mechanism(g)
p <- enumerate_data(g); pt <- study_task(p,g)$task
nu <- population_truth_functions(p,g);theta <- exact_truth(g)
stopifnot(max(abs(theta-backward_truth(g)))<1e-12)
audits <- equations <- list()
for(time in 1:2) {
 H <- pt$vars$history('A',time);A <- pt$vars$A[[time]]
 hkey <- cell_key(pt$natural[,H,drop=FALSE])
 bkey <- cell_key(pt$natural[,c(H,A),drop=FALSE])
 ckey <- cell_key(p[,c(paste0('C',time+1L,'_covariate'),paste0('Y',time+1L)),drop=FALSE])
 m <- p[[paste0('R',time+1L)]]==1
 vkey <- paste(hkey[m],ckey[m],sep='|')
 w <- p$probability;y <- ifelse(m,p[[paste0('Y',time+1L)]],0)
 omega <- apply(nu$ratios[,seq_len(time),drop=FALSE],1,prod)
 previous <- if(time==1L) rep(1,nrow(p)) else nu$ratios[,1]
 residual <- rep(-1,nrow(p));residual[m]<-nu$beta[m,time]-1
 br <- conditional_mean(residual,w,bkey)
 ar <- conditional_mean(nu$lambda[m,time]-omega[m]*y[m],w[m],vkey)
 bridge_mean <- sum(w*omega*y*nu$beta[,time]);adjoint_mean<-sum(w*nu$lambda[,time])
 equations[[time]]<-data.table(time=time,outcome_time=time+1L,theta=theta[time],
  forward_backward_error=theta[time]-backward_truth(g)[time],
  true_bridge_equation_max_error=max(abs(br)),true_adjoint_equation_max_error=max(abs(ar)),
  bridge_representation=bridge_mean,adjoint_representation=adjoint_mean,
  bridge_representation_error=bridge_mean-theta[time],adjoint_representation_error=adjoint_mean-theta[time])
 for(history in unique(hkey)) {
  rows<-which(hkey==history);observed<-rows[m[rows]]
  bs<-sort(unique(bkey[rows]));cs<-sort(unique(ckey[observed]))
  joint<-as.matrix(Matrix::sparseMatrix(i=match(bkey[observed],bs),j=match(ckey[observed],cs),
   x=w[observed],dims=c(length(bs),length(cs))))
  marginal<-vapply(bs,function(cell)sum(w[rows[bkey[rows]==cell]]),0)
  T<-joint/marginal;Ts<-t(joint)/colSums(joint)
  d<-svd(T)$d;dstar<-svd(Ts)$d
  rank<-sum(d>1e-10*max(d));rankstar<-sum(dstar>1e-10*max(dstar))
  audits[[length(audits)+1L]]<-data.table(time=time,history=history,
   natural_history_probability=sum(w[rows]),policy_history_probability=sum(w[rows]*previous[rows]),
   current_treatment_values=length(bs),complete_health_values=length(cs),
   bridge_rank=rank,bridge_nullity=length(cs)-rank,
   adjoint_rank=rankstar,adjoint_nullity=length(bs)-rankstar,
   smallest_bridge_singular_value=min(d),bridge_condition_number=max(d)/min(d))
 }
 if(time==2L) {
  # Demonstrate two distinct valid adjoints on a positive-policy-probability history.
  history<-unique(hkey[previous>0])[1];rows<-which(hkey==history);observed<-rows[m[rows]]
  bs<-sort(unique(bkey[rows]));cs<-sort(unique(ckey[observed]))
  joint<-as.matrix(Matrix::sparseMatrix(i=match(bkey[observed],bs),j=match(ckey[observed],cs),
   x=w[observed],dims=c(length(bs),length(cs))))
  Ts<-t(joint)/colSums(joint)
  decomposition<-eigen(crossprod(Ts),symmetric=TRUE)
  difference<-decomposition$vectors[,ncol(decomposition$vectors)]
  difference<-difference/max(abs(difference))
  alt<-nu$lambda[,time];alt[rows]<-alt[rows]+difference[match(bkey[rows],bs)]
  alt_error<-conditional_mean(alt[m]-omega[m]*y[m],w[m],vkey)
  saved<-readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
  estimated_beta<-saved$fitted_functions[['sdr/beta1_lambda1_ratio1_m1/2']]$beta[,2]
  rr<-rep(-1,nrow(p));rr[m]<-estimated_beta[m]-1
  alternative<-data.table(time=2L,maximum_pointwise_adjoint_difference=max(abs(alt-nu$lambda[,time])),
   conditional_equation_max_error=max(abs(alt_error)),
   parameter_mean_difference=sum(w*(alt-nu$lambda[,time])),
   correction_mean_difference=-sum(w*(alt-nu$lambda[,time])*rr))
  stopifnot(max(abs(alt_error))<1e-10,abs(alternative$parameter_mean_difference)<1e-10,
   abs(alternative$correction_mean_difference)<1e-10)
  fwrite(alternative,file.path(out,'alternative-adjoint-check.csv'))
 }
}
a<-rbindlist(audits);eq<-rbindlist(equations)
stopifnot(all(a$bridge_nullity==0L),all(a$adjoint_nullity==4L),
 all(eq$true_bridge_equation_max_error<1e-10),all(eq$true_adjoint_equation_max_error<1e-10),
 all(abs(eq$bridge_representation_error)<1e-12),all(abs(eq$adjoint_representation_error)<1e-12))
fwrite(a,file.path(out,'full-history-operator-audit.csv'))
fwrite(eq,file.path(out,'true-equation-and-parameter-audit.csv'))
constraints<-data.table(condition=c('R1=R2=1','coded deterministic treatment at R_t=0',
 'policy support on reachable histories','positive final measurement',
 'bridge exists on all reachable histories','adjoint exists on all reachable histories'),
 pass=c(all(p$R1==1 & p$R2==1),all(g$g[,,1,1]==1)&all(g$g[,,1,-1]==0),
 all(g$gd[g$g==0]==0)&all(a$natural_history_probability[a$policy_history_probability>0]>0),
 min(g$measurement[,2,,])>0,all(a$bridge_nullity==0),all(a$adjoint_rank==a$complete_health_values)))
stopifnot(all(constraints$pass))
fwrite(constraints,file.path(out,'assumption-audit.csv'))
print(eq,digits=15)
print(a[,.(histories=.N,bridge_rank_min=min(bridge_rank),bridge_rank_max=max(bridge_rank),
 bridge_nullity=max(bridge_nullity),adjoint_nullity=max(adjoint_nullity),
 min_singular_value=min(smallest_bridge_singular_value),max_condition_number=max(bridge_condition_number)),by=time])
print(alternative,digits=15)
cat('Final measurement probability minimum:',min(g$measurement[,2,,]),'\n')
cat('All reachable full histories pass. Non-visited treatment histories have probability zero under this DGP and policy.\n')
cat('Assumptions 2 and 4 follow from the independent structural draws and absence of treatment from the measurement equation.\n')
