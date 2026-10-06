# Population-only check using the complete observed history, independently
# rebuilding the conditional bridge operator from enumerated probabilities.
source(file.path(Sys.getenv('STUDY_SOURCE','simulations/observed_history_corrected'),'study.R'))
root <- Sys.getenv('SIM_OUTPUT')
blocks <- examples <- list()
for (mechanism in c('binary_longitudinal','discrete_dose')) {
  g <- make_mechanism(mechanism,dose_max=3L,visits='missing_both')
  p <- enumerate_data(g);pt <- study_task(p,g)$task
  nu <- population_truth_functions(p,g)
  for (t in 1:2) {
    H <- pt$vars$history('A',t);A <- pt$vars$A[[t]]
    health <- c(paste0('C',t+1L,'_covariate'),paste0('Y',t+1L))
    measured <- p[[paste0('R',t+1L)]]==1
    data <- data.table(H=cell_key(pt$natural[,H,drop=FALSE]),
      A=cell_key(pt$natural[,A,drop=FALSE]),C=cell_key(p[,health,drop=FALSE]),
      measured,current_R=p[[paste0('R',t)]],probability=p$probability,
      beta=nu$beta[,t],baseline=p$C1_baseline)
    # All upstream observed history is retained in H; no generating-variable
    # reduction and no sample-based rank calculation is used.
    for (history in unique(data$H)) {
      d <- data[H==history]
      action <- sort(unique(d$A));c_values <- sort(unique(d$C[d$measured]))
      denominator <- d[,.(mass=sum(probability)),by=A]
      numerator <- d[measured==TRUE,.(mass=sum(probability)),by=.(A,C)]
      T <- matrix(0,length(action),length(c_values),dimnames=list(action,c_values))
      for (j in seq_len(nrow(numerator))) T[numerator$A[j],numerator$C[j]] <-
        numerator$mass[j]/denominator$mass[match(numerator$A[j],denominator$A)]
      spectrum <- svd(T)$d
      rank <- sum(spectrum>max(spectrum)*1e-12)
      beta <- vapply(c_values,function(c) unique(d$beta[d$measured & d$C==c])[1L],numeric(1))
      stopifnot(max(abs(T%*%beta-1))<1e-12)
      blocks[[length(blocks)+1L]] <- data.table(mechanism,time=t,H=history,current_R=d$current_R[1L],
        history_probability=sum(d$probability),treatment_values=nrow(T),health_values=ncol(T),rank,
        unique_bridge=rank==ncol(T),smallest_reported_singular_value=min(spectrum),
        largest_singular_value=max(spectrum),selected_solution_error=max(abs(T%*%beta-1)))
      if (t==2 && d$current_R[1L]==0) {
        # Compare the selected constant solution with the inverse missingness
        # probability. Both satisfy the equation and have finite inverse-expit
        # coefficients, proving nonuniqueness even in that parameterization.
        q <- d[d$measured]
        index <- which(data$H==history & data$measured)
        original_index <- which(data$H==history & data$measured)
        state <- state_index(p,t)
        ii <- which(data$H==history & measured)
        health_index <- 1L+p[[health[1]]][ii]+2L*p[[health[2]]][ii]
        inverse <- 1/g$measurement[cbind(p$C1_baseline[ii]+1L,t,state[ii],health_index)]
        inverse_values <- vapply(c_values,function(c) inverse[match(c,data$C[ii])],numeric(1))
        stopifnot(all(beta>1),all(inverse_values>1),max(abs(T%*%inverse_values-1))<1e-12,
          max(abs(beta-inverse_values))>1e-3)
        if (sum(vapply(examples,function(x)x$mechanism[1L]==mechanism,logical(1)))==0L)
          examples[[length(examples)+1L]] <- data.table(mechanism,H=history,baseline=d$baseline[1L],
            C=c_values,conditional_coefficient=as.numeric(T),selected_constant=beta,
            inverse_probability_bridge=inverse_values,
            selected_link_coefficient=-log(beta-1),inverse_link_coefficient=-log(inverse_values-1),
            selected_equation=sum(T*beta),inverse_equation=sum(T*inverse_values))
      }
    }
  }
}
blocks <- rbindlist(blocks);examples<-rbindlist(examples)
summary <- blocks[,.(histories=.N,treatment_values=paste(sort(unique(treatment_values)),collapse=','),
 health_values=paste(sort(unique(health_values)),collapse=','),rank=paste(sort(unique(rank)),collapse=','),
 unique_at_every_history=all(unique_bridge),minimum_singular_value=min(smallest_reported_singular_value)),
 by=.(mechanism,time,current_R)]
stopifnot(all(blocks$unique_bridge[blocks$current_R==1]),!any(blocks$unique_bridge[blocks$current_R==0]))
fwrite(blocks,file.path(root,'bridge-uniqueness-full-history.csv'))
fwrite(summary,file.path(root,'bridge-uniqueness-summary.csv'))
fwrite(examples,file.path(root,'two-positive-bridge-solutions.csv'))
print(summary);print(examples)
cat('All full observed histories checked. Bridges unique at visited treatment times; not unique at missed time-2 visits, even within inverse-expit.\n')
