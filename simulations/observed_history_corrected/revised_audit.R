# Deterministic calculations for the revised two-mechanism plan only.
# No samples are drawn and no estimator or learner is fitted.
source('simulations/observed_history_corrected/study.R')
out <- 'results/observed-history-corrected/revised-design'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
audits <- truths <- list()
for (name in c('binary_longitudinal', 'discrete_dose')) {
  g <- make_mechanism(name, visits = 'scheduled', dose_max = 3L)
  audits[[name]] <- require_valid_mechanism(g)
  stopifnot(all(audits[[name]]$operator_rank == nrow(g$C)))
  forward <- exact_truth(g); backward <- backward_truth(g)
  p <- enumerate_data(g); nu <- population_truth_functions(p, g)
  bridge <- adjoint <- numeric(g$tau)
  omega <- rep(1, nrow(p))
  for (t in seq_len(g$tau)) {
    omega <- omega * nu$ratios[, t]
    measured <- p[[paste0('R', t + 1L)]] == 1
    bridge[t] <- sum(p$probability[measured] * omega[measured] *
      p[[paste0('Y', t + 1L)]][measured] * nu$beta[measured, t])
    adjoint[t] <- sum(p$probability * nu$lambda[, t])
  }
  stopifnot(max(abs(forward - backward)) < 1e-12,
    max(abs(forward - bridge)) < 1e-12, max(abs(forward - adjoint)) < 1e-12)
  truths[[name]] <- data.frame(mechanism = name, outcome_time = seq_len(g$tau) + 1L,
    dose_max = g$dose_max, truth = forward, backward_truth = backward,
    bridge_representation = bridge, adjoint_representation = adjoint)
}
write.csv(do.call(rbind, audits), file.path(out, 'assumption_audit.csv'), row.names = FALSE)
write.csv(do.call(rbind, truths), file.path(out, 'truth.csv'), row.names = FALSE)
print(do.call(rbind, truths), digits = 15)
cat('Both revised mechanisms pass; bridge operators have full column rank.\n')
cat('Independent forward/backward sums and both Theorem 2 representations agree within 1e-12.\n')
cat('No estimator simulations or generation pilots were run.\n')
