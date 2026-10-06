# Audit the saved unpenalized bridge against empirical and population equations.
# No new data, learner fit, or known-function performance simulation is run.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
inputs <- readRDS(file.path(out, 'scoring-inputs.rds'))
fit <- readRDS(file.path(out, 'zero-penalty-bridge-fit.rds'))
points <- fread(file.path(out, 'bridge-predictions.csv'))[source == 'large check']
saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
p <- saved$population
g <- make_mechanism('discrete_dose', dose_max = 3L)
pt <- study_task(p, g)$task
nu <- population_truth_functions(p, g)
H <- pt$vars$history('A', 2L); treatment <- pt$vars$A[[2L]]
Bpop <- as.matrix(pt$natural[, c(H, treatment), drop = FALSE])
Vpop <- as.matrix(cbind(pt$natural[, H, drop = FALSE],
                       p[, c('C3_covariate', 'Y3'), drop = FALSE]))
stopifnot(identical(colnames(inputs$B), colnames(Bpop)))
zkey <- cell_key(inputs$B); zlevels <- sort(unique(zkey))
xlevels <- fit$target_levels
zi <- match(zkey, zlevels)
xi <- match(inputs$beta_key[inputs$observed], xlevels)
counts <- tabulate(zi, length(zlevels))
joint <- as.matrix(Matrix::sparseMatrix(i = zi[inputs$observed], j = xi, x = 1,
  dims = c(length(zlevels), length(xlevels))))
empirical_A <- joint / counts
weights <- counts / sum(counts)
pop_zi <- match(cell_key(Bpop), zlevels)
pop_xi <- match(cell_key(Vpop), xlevels)
observed <- p$R3 == 1L
stopifnot(!anyNA(pop_zi), !anyNA(pop_xi[observed]), !anyNA(xi))
pop_bprob <- as.numeric(Matrix::sparseMatrix(i = pop_zi, j = rep(1L, nrow(p)),
  x = p$probability, dims = c(length(zlevels), 1L)))
pop_joint <- as.matrix(Matrix::sparseMatrix(i = pop_zi[observed], j = pop_xi[observed],
  x = p$probability[observed], dims = c(length(zlevels), length(xlevels))))
population_A <- pop_joint / pop_bprob
# The saved distribution represents observed data: unmeasured C3 is NA.
# Independently obtain the bridge from the generating measurement equation.
first <- match(xlevels, cell_key(Vpop))
health_index <- 1L + p$C3_covariate[first] + 2L * p$Y3[first]
true_beta <- 1 / g$measurement[cbind(p$C1_baseline[first] + 1L, 2L,
  state_index(p[first, ], 2L), health_index)]
stopifnot(max(abs(true_beta[pop_xi[observed]] - nu$beta[observed, 2L])) < 1e-12,
          all(true_beta >= fit$tuning$lower & true_beta <= fit$tuning$upper))

# Exact conditional sampling variance of the bridge equation at the truth.
# Every treatment/history row is a multinomial sample of four measured
# health configurations plus an unmeasured category.
equation_variance <- pmax(0, as.vector(population_A %*% true_beta^2) - 1) / counts
equation_se <- sqrt(equation_variance)
empirical_true_lhs <- as.vector(empirical_A %*% true_beta)
standardized_equation_error <- (empirical_true_lhs - 1) / equation_se
# A first-order sampling-SD approximation for the unpenalized coefficients.
# The probabilities enter only this diagnostic, never either saved fit.
weighted_population_A <- population_A * sqrt(weights)
linear_map <- qr.coef(qr(weighted_population_A, tol = 1e-12), diag(sqrt(weights)))
beta_sampling_sd <- sqrt(rowSums((linear_map * rep(equation_se, each = length(xlevels)))^2))
population_solution <- qr.coef(qr(population_A, tol = 1e-12), rep(1, nrow(population_A)))
stopifnot(max(abs(population_solution - true_beta)) < 1e-10)

selected <- points$estimate[match(xlevels, points$cell)]
raw <- fit$coefficients[1L] + fit$coefficients[-1L]
direct <- fit$direct_coefficients
bounded <- pmax(fit$tuning$lower, pmin(fit$tuning$upper, raw))
scenarios <- list('True bridge (evaluation only)' = true_beta,
  'CV-selected positive penalty' = selected,
  'Zero L1 penalty, original bounds' = bounded,
  'Zero L1 penalty, no prediction bounds' = raw,
  'Independent unpenalized least squares' = direct)
summary <- rbindlist(lapply(names(scenarios), function(name) {
  beta <- scenarios[[name]]
  empirical <- as.vector(empirical_A %*% beta) - 1
  population <- as.vector(population_A %*% beta) - 1
  data.table(scenario = name,
    empirical_equation_rmse = sqrt(sum(weights * empirical^2)),
    empirical_equation_maximum_absolute_residual = max(abs(empirical)),
    population_equation_rmse = sqrt(sum(pop_bprob * population^2)),
    population_equation_maximum_absolute_residual = max(abs(population)))
}))

# The independent LS estimate minus the truth must equal the LS projection
# of the empirical equation's sampling discrepancy at the true bridge.
# This is an algebraic decomposition, not a learner trained on true values.
weighted_A <- empirical_A * sqrt(weights)
discrepancy <- 1 - as.vector(empirical_A %*% true_beta)
projected_discrepancy <- qr.coef(qr(weighted_A, tol = 1e-12), sqrt(weights) * discrepancy)
identity_error <- max(abs((direct - true_beta) - projected_discrepancy))
span <- data.table(reachable_population_beta_values = length(true_beta),
  represented_training_beta_values = length(xlevels),
  empirical_operator_rank = qr(weighted_A)$rank,
  true_minimum = min(true_beta), true_maximum = max(true_beta),
  true_population_equation_maximum_residual = max(abs(population_A %*% true_beta - 1)),
  maximum_population_solution_difference_from_generating_truth = max(abs(population_solution - true_beta)),
  maximum_sampling_discrepancy_decomposition_error = identity_error,
  first_order_predicted_probability_weighted_beta_rmse = sqrt(weighted.mean(beta_sampling_sd^2, colSums(pop_joint))),
  mean_squared_standardized_empirical_equation_error_at_truth = mean(standardized_equation_error^2),
  maximum_absolute_standardized_empirical_equation_error_at_truth = max(abs(standardized_equation_error)),
  treatment_history_groups_with_absolute_standardized_error_over_3 = sum(abs(standardized_equation_error) > 3))
stopifnot(identity_error < 1e-9,
          span$true_population_equation_maximum_residual < 1e-12)

point_table <- as.data.table(p[first, c('C1_baseline', 'A1_policy', 'A1_other',
                                      'C2_covariate', 'Y2', 'C3_covariate', 'Y3')])
point_table[, `:=`(cell = xlevels, history = cell_key(Vpop[first, H, drop = FALSE]),
  true = true_beta, selected_estimate = selected, zero_raw_estimate = raw,
  zero_bounded_estimate = bounded, independent_zero_estimate = direct,
  zero_raw_error = raw - true_beta, zero_bounded_error = bounded - true_beta,
  training_measured = colSums(joint), expected_training_measured = sum(counts) * colSums(pop_joint),
  measured_population_probability = colSums(pop_joint) / sum(pop_joint),
  first_order_sampling_sd = beta_sampling_sd,
  unpenalized_error_divided_by_first_order_sd = (direct - true_beta) / beta_sampling_sd)]
point_table <- point_table[order(-abs(zero_raw_error))]
worst_history <- point_table$history[1L]
pop_history <- cell_key(Vpop[, H, drop = FALSE])
pop_zfirst <- match(zlevels, cell_key(Bpop))
row_history <- cell_key(Bpop[pop_zfirst, H, drop = FALSE])
target_history <- pop_history[first]
ri <- which(row_history == worst_history); ci <- which(target_history == worst_history)
worst <- list()
for (row in ri) {
  for (column in ci) {
    decoded <- as.data.table(p[pop_zfirst[row], c('A2_policy', 'A2_other')])
    decoded[, `:=`(history = worst_history, cell = xlevels[column],
      C3_covariate = p$C3_covariate[first[column]], Y3 = p$Y3[first[column]],
      training_treatment_history_n = counts[row], measured_health_n = joint[row, column],
      empirical_joint_probability = empirical_A[row, column],
      population_joint_probability = population_A[row, column],
      true_bridge = true_beta[column], unpenalized_bridge = direct[column],
      true_bridge_empirical_equation_lhs = as.vector(empirical_A[row, , drop = FALSE] %*% true_beta),
      true_bridge_population_equation_lhs = as.vector(population_A[row, , drop = FALSE] %*% true_beta),
      equation_sampling_sd_at_truth = equation_se[row],
      standardized_empirical_equation_error_at_truth = standardized_equation_error[row],
      unpenalized_bridge_empirical_equation_lhs = as.vector(empirical_A[row, , drop = FALSE] %*% direct),
      unpenalized_bridge_population_equation_lhs = as.vector(population_A[row, , drop = FALSE] %*% direct))]
    worst[[length(worst) + 1L]] <- decoded
  }
}
fwrite(summary, file.path(out, 'zero-penalty-equation-comparison.csv'))
fwrite(span, file.path(out, 'zero-penalty-span-and-decomposition.csv'))
fwrite(point_table, file.path(out, 'zero-penalty-point-diagnosis.csv'))
fwrite(rbindlist(worst), file.path(out, 'zero-penalty-worst-history-equations.csv'))
print(summary, digits = 9)
print(span, digits = 9)
print(point_table[1:8, .(true, selected_estimate, zero_raw_estimate,
  zero_raw_error, training_measured, expected_training_measured,
  measured_population_probability, first_order_sampling_sd,
  unpenalized_error_divided_by_first_order_sd)], digits = 9)
cat('Worst-history total training people:', sum(counts[ri]), '\n')
cat('Worst-history treatment-group training counts:', paste(counts[ri], collapse = ', '), '\n')
cat('Worst-history measured health counts:', paste(colSums(joint[ri, ci]), collapse = ', '), '\n')
cat('The true function is evaluated, not fitted. All bridge estimates are the existing empirical fits.\n')
