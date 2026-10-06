# Empirical conditional equations on the original n=4,000 sample.
# Uses saved fits from outside each observation's estimator fold; no refitting.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
g <- make_mechanism('discrete_dose', dose_max = 3L)
set.seed(saved$design$seed)
d <- draw_data(saved$design$n, g)
p <- saved$population
columns <- names(d)
pkey <- cell_key(p[, columns, drop = FALSE])
index <- match(cell_key(d), pkey)
stopifnot(!anyNA(index))
beta <- lambda <- omega <- rep(NA_real_, nrow(d))
fold <- integer(nrow(d))
for (j in seq_along(saved$folds)) {
  fitted <- saved$fitted_functions[[paste('sdr', 'beta1_lambda1_ratio1_m1', j, sep = '/')]]
  # Repeated observed rows in the enumerated distribution differ only in
  # unobserved complete health; the saved observable predictions must agree.
  check <- data.table(cell = pkey, beta = fitted$beta[, 2L],
    lambda = fitted$lambda[, 2L], omega = apply(fitted$ratios, 1L, prod))[,
    .(beta_range = diff(range(beta)), lambda_range = diff(range(lambda)),
      omega_range = diff(range(omega))), by = cell]
  stopifnot(all(check$beta_range == 0), all(check$lambda_range == 0), all(check$omega_range == 0))
  rows <- saved$folds[[j]]$validation_set
  beta[rows] <- fitted$beta[index[rows], 2L]
  lambda[rows] <- fitted$lambda[index[rows], 2L]
  omega[rows] <- apply(fitted$ratios[index[rows], , drop = FALSE], 1L, prod)
  fold[rows] <- j
}
stopifnot(all(fold > 0L), all(is.finite(beta)), all(is.finite(lambda)), all(is.finite(omega)))
d$fold <- fold
d$beta <- beta
d$lambda <- lambda
d$omega <- omega
d$bridge_lhs <- ifelse(d$R3 == 1L, beta, 0)
d$adjoint_response <- ifelse(d$R3 == 1L, omega * d$Y3, 0)
fwrite(d, file.path(out, 'sample-equation-observations.csv'))

# Constants and missingness indicators add no conditioning information in
# this DGP (R1=R2=1, baseline covariate and outcome are zero).
H <- c('C1_baseline', 'A1_policy', 'A1_other', 'C2_covariate', 'Y2')
B <- c(H, 'A2_policy', 'A2_other')
V <- c(H, 'C3_covariate', 'Y3')
summaries <- bridges <- adjoints <- list()
for (j in c(0L, 1:3)) {
  scope <- if (j == 0L) 'All 4,000 observations, cross-fitted' else paste0('Training sample ', j, ': held-out observations')
  rows <- if (j == 0L) seq_len(nrow(d)) else which(d$fold == j)
  current <- as.data.table(d[rows, ])
  b <- current[, .(observations = .N, measured = sum(R3), lhs = mean(bridge_lhs), rhs = 1), by = B]
  b[, `:=`(residual = lhs - rhs, scope = scope, fold = j)]
  a <- current[R3 == 1L, .(observations = .N, lhs = mean(lambda), rhs = mean(adjoint_response)), by = V]
  a[, `:=`(residual = lhs - rhs, scope = scope, fold = j)]
  bridges[[length(bridges) + 1L]] <- b
  adjoints[[length(adjoints) + 1L]] <- a
  for (kind in c('bridge', 'adjoint')) {
    cells <- if (kind == 'bridge') b else a
    summaries[[length(summaries) + 1L]] <- data.table(
      n = nrow(d), seed = saved$design$seed, time = 2L, fold = j, scope = scope,
      equation = kind, observations = sum(cells$observations), groups = nrow(cells),
      mean_lhs = weighted.mean(cells$lhs, cells$observations),
      mean_rhs = weighted.mean(cells$rhs, cells$observations),
      mean_residual = weighted.mean(cells$residual, cells$observations),
      conditional_equation_rmse = sqrt(weighted.mean(cells$residual^2, cells$observations)),
      maximum_absolute_conditional_residual = max(abs(cells$residual)),
      minimum_group_size = min(cells$observations), median_group_size = median(cells$observations),
      groups_with_one_observation = sum(cells$observations == 1L),
      minimum_lhs = min(cells$lhs), maximum_lhs = max(cells$lhs))
  }
}
fwrite(rbindlist(bridges), file.path(out, 'sample-bridge-equation-values.csv'))
fwrite(rbindlist(adjoints), file.path(out, 'sample-adjoint-equation-values.csv'))
summary <- rbindlist(summaries)
fwrite(summary, file.path(out, 'sample-equation-summary.csv'))
print(summary[, .(fold, equation, observations, groups, mean_lhs, mean_rhs,
  mean_residual, conditional_equation_rmse, maximum_absolute_conditional_residual)], digits = 8)
cat('The adjoint right side uses estimated, out-of-fold treatment ratios.\n')
cat('These are empirical conditional means on evaluation rows, not population-DGP expectations.\n')
cat('At time 1 R2=1 and beta1=1, so every empirical bridge left side is exactly one.\n')
cat('No adjoint is fitted at time 1 in the estimator.\n')
