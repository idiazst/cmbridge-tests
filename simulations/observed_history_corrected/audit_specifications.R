source("simulations/observed_history_corrected/study.R")
answer <- list()
for (name in c("binary_point", "binary_longitudinal", "discrete_dose")) {
  g <- make_mechanism(name); p <- enumerate_data(g); pt <- study_task(p, g)$task
  nu <- population_truth_functions(p, g); s <- g$tau; w <- p$probability
  measured <- p[[paste0("R", s + 1L)]] == 1
  H <- pt$vars$history("A", s); A <- pt$vars$A[[s]]
  B <- cell_key(pt$natural[, c(H, A), drop = FALSE])
  V <- cell_key(cbind(pt$natural[measured, H, drop = FALSE],
    p[measured, c(paste0("C", s + 1L, "_covariate"), paste0("Y", s + 1L)), drop = FALSE]))
  probability <- conditional_mean(as.numeric(measured), w, B)
  omega <- apply(nu$ratios[, seq_len(s), drop = FALSE], 1, prod)
  loading <- conditional_mean(omega[measured] * p[[paste0("Y", s + 1L)]][measured], w[measured], V)
  terminal <- numeric(nrow(p)); terminal[measured] <- nu$beta[measured, s] * p[[paste0("Y", s + 1L)]][measured]
  m <- conditional_mean(terminal, w, B)
  variations <- vapply(seq_len(s), function(t) {
    key <- cell_key(pt$natural[, pt$vars$A[[t]], drop = FALSE])
    restricted <- conditional_mean(nu$ratios[, t], w, key)
    sum(w * (nu$ratios[, t] - unname(restricted[key]))^2)
  }, 0)
  row <- data.table(mechanism = name, main_horizon = s,
    measurement_given_A_H_range = diff(range(probability)),
    adjoint_target_given_H_C_range = diff(range(loading)),
    terminal_regression_range = diff(range(m)),
    ratio_history_omission_mse = sum(variations))
  stopifnot(row$measurement_given_A_H_range > 1e-8,
    row$adjoint_target_given_H_C_range > 1e-8, row$terminal_regression_range > 1e-8,
    row$ratio_history_omission_mse > 1e-8)
  answer[[name]] <- row
}
fwrite(rbindlist(answer), "results/observed-history-corrected/specification_audit.csv")
print(rbindlist(answer))
cat("Omitted variation verified by deterministic population calculations; no estimator simulation with known functions.\n")
