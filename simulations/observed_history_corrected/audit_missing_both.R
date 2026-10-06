# Exact checks of the revised DGP. No samples or fitted estimators are produced.
source('simulations/observed_history_corrected/study.R')
out <- 'results/observed-history-corrected/missing-both-design'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
blocks <- truths <- missing <- equations <- classes <- alternatives <- list()
for (name in c('binary_longitudinal','discrete_dose')) {
  g <- make_mechanism(name, visits='missing_both', dose_max=3L)
  aa <- require_valid_mechanism(g)
  stopifnot(all(aa$operator_rank[!aa$no_visit] == nrow(g$C)),
            all(aa$operator_rank[aa$no_visit] == 1L),
            all(aa$minimum_singular_value > 0))
  blocks[[name]] <- aa
  forward <- exact_truth(g); backward <- backward_truth(g)
  p <- enumerate_data(g)
  pt <- study_task(p,g)$task
  nu <- population_truth_functions(p,g)
  omega <- rep(1,nrow(p)); bridge <- adjoint <- numeric(g$tau)
  alternative_beta <- nu$beta
  for(t in seq_len(g$tau)) {
    omega <- omega*nu$ratios[,t]
    measured <- p[[paste0('R',t+1L)]]==1L
    H <- pt$vars$history('A',t); A <- pt$vars$A[[t]]
    B <- cell_key(pt$natural[,c(H,A),drop=FALSE])
    V <- cell_key(cbind(pt$natural[,H,drop=FALSE],p[,c(paste0('C',t+1L,'_covariate'),paste0('Y',t+1L)),drop=FALSE]))
    bridge_residual <- rep(-1,nrow(p))
    bridge_residual[measured] <- nu$beta[measured,t]-1
    bm <- conditional_mean(bridge_residual,p$probability,B)
    am <- conditional_mean(nu$lambda[measured,t]-omega[measured]*p[[paste0('Y',t+1L)]][measured],p$probability[measured],V[measured])
    stopifnot(max(abs(bm))<1e-10,max(abs(am))<1e-10)
    equations[[length(equations)+1L]] <- data.table(mechanism=name,time=t,
      bridge_histories=length(bm),adjoint_health_histories=length(am),
      bridge_max_error=max(abs(bm)),adjoint_max_error=max(abs(am)))
    for(kind in c('bridge','adjoint','ratio')) {
      key <- if(kind=='bridge')V[measured] else B
      value <- if(kind=='bridge')nu$beta[measured,t] else
        if(kind=='adjoint')nu$lambda[,t] else nu$ratios[,t]
      span <- data.table(cell=key,value=value)[,.(spread=max(value)-min(value)),by=cell]
      stopifnot(max(span$spread)<1e-10)
      classes[[length(classes)+1L]] <- data.table(mechanism=name,time=t,function_name=kind,
        joint_combinations=nrow(span),maximum_same_input_spread=max(span$spread),
        full_joint_category_class_contains_solution=TRUE)
    }
    bridge[t] <- sum(p$probability[measured]*omega[measured]*p[[paste0('Y',t+1L)]][measured]*nu$beta[measured,t])
    adjoint[t] <- sum(p$probability*nu$lambda[,t])
    # The inverse measurement probability is another valid bridge on the
    # missed-visit histories, distinct from the constant selected solution.
    if(t==2L) {
      ii <- which(measured & p$R2==0L)
      h <- state_index(p,t); c <- 1L+p$C3_covariate[ii]+2L*p$Y3[ii]
      alternative_beta[ii,t] <- 1/g$measurement[cbind(p$C1_baseline[ii]+1L,t,h[ii],c)]
      rr <- rep(-1,nrow(p));rr[measured] <- alternative_beta[measured,t]-1
      check <- conditional_mean(rr,p$probability,B)
      alternative_truth <- sum(p$probability[measured]*omega[measured]*p$Y3[measured]*alternative_beta[measured,t])
      stopifnot(max(abs(check))<1e-10,abs(alternative_truth-forward[t])<1e-12,
        all(nu$beta[measured,t]>=1 & nu$beta[measured,t]<=6),
        max(abs(alternative_beta[ii,t]-nu$beta[ii,t]))>1e-3)
      alternatives[[name]] <- data.table(mechanism=name,
        maximum_solution_difference=max(abs(alternative_beta[ii,t]-nu$beta[ii,t])),
        alternative_bridge_equation_error=max(abs(check)),
        parameter_difference=alternative_truth-forward[t])
    }
  }
  stopifnot(max(abs(forward-backward))<1e-12,max(abs(forward-bridge))<1e-12,
            max(abs(forward-adjoint))<1e-12,all(g$beta>=1 & g$beta<=6))
  truths[[name]] <- data.table(mechanism=name,outcome_time=2:3,truth=forward,
    backward_truth=backward,bridge_representation=bridge,adjoint_representation=adjoint)
  missing[[name]] <- data.table(mechanism=name,time=1:3,
    natural_missing_percent=100*apply(g$natural_states,2,function(z)sum(z[,1L])),
    policy_missing_percent=100*apply(g$policy_states,2,function(z)sum(z[,1L])))
  # Deterministic no-visit treatment and the baseline-outcome transition
  # are checked on the entire finite observed support, including policy
  # histories (which are covered by natural support).
  stopifnot(all(p$A2_policy[p$R2==0]==0),all(p$A2_other[p$R2==0]==0),
    all(p$Y3[p$R2==0 & p$R3==1]==p$Y1[p$R2==0 & p$R3==1]),
    sum(p$probability[p$R2==0 & p$R3==1])>0,
    sum(p$probability[p$R2==1 & p$R3==0])>0,
    sum(p$probability[p$R2==0 & p$R3==0])>0)
}
fwrite(rbindlist(blocks),file.path(out,'assumption_audit.csv'))
fwrite(rbindlist(truths),file.path(out,'truth.csv'))
fwrite(rbindlist(missing),file.path(out,'missingness.csv'))
fwrite(rbindlist(equations),file.path(out,'full-history-equations.csv'))
fwrite(rbindlist(classes),file.path(out,'function-classes.csv'))
fwrite(rbindlist(alternatives),file.path(out,'alternative-bridge-solutions.csv'))
print(rbindlist(truths),digits=15);print(rbindlist(missing),digits=10)
print(rbindlist(equations));print(rbindlist(classes));print(rbindlist(alternatives))
cat('Both revised mechanisms satisfy the checked paper conditions.\n')
cat('The bridge operator has rank four at a measured treatment visit, and rank one at a missed visit.\n')
cat('Forward/backward truths and both identifying representations agree within 1e-12.\n')
cat('No samples were generated and no estimators were fitted.\n')
