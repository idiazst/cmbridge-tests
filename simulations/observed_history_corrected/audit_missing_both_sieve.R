# Check saved sieve models without fitting another learner or generating data.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
root <- Sys.getenv('SIM_OUTPUT')
mechanism <- Sys.getenv('DIAGNOSTIC_MECHANISM', 'discrete_dose')
replicate <- as.integer(Sys.getenv('DIAGNOSTIC_REPLICATE', '6'))
tag <- if (mechanism == 'binary_longitudinal') 'binary-dose' else 'discrete-dose'
out <- file.path(root, sprintf('diagnostics/%s-n4000-r%03d-sieve', tag, replicate))
saved <- readRDS(file.path(out, 'recovered-sieve-fits.rds'))
original <- readRDS(file.path(root, sprintf('jobs/%s-n4000-r%03d.rds', mechanism, replicate)))
d <- saved$data
g <- make_mechanism(mechanism, dose_max = max(2L, saved$design$dose_max), visits = saved$design$visits)
prepared <- study_task(d, g, folds = saved$folds, learner_folds = saved$learner_folds,
  learner_groups = saved$design$learner_groups, balance_measurement = saved$design$balance_measurement)
task <- prepared$task
H <- task$vars$history('A', 2L); A <- task$vars$A[[2L]]
V <- as.matrix(cbind(task$natural[, H, drop = FALSE], d[, c('C3_covariate', 'Y3')]))
B <- as.matrix(task$natural[, c(H, A), drop = FALSE])
Hkeys <- cell_key(task$natural[, H, drop = FALSE])
points <- fread(file.path(out, 'bridge-predictions.csv'))
summary <- partitions <- histories <- direct_checks <- list()
for (j in 1:3) {
  model <- saved$models[[j]]
  rows <- saved$folds[[j]]$training_set
  M <- d$R3[rows]; n <- length(rows)
  xkeys <- cell_key(V[rows[M == 1L], , drop = FALSE])
  zkeys <- cell_key(B[rows, , drop = FALSE])
  xi <- match(xkeys, model$target_spec$levels)
  zi <- match(zkeys, model$instrument_spec$levels)
  # Independently construct the individual-row indicator designs. Use a
  # Cholesky whitening transform and augmented QR least squares, rather
  # than the package's grouped normal equations and matrix inversion.
  X <- matrix(0, n, length(model$coefficients))
  X[M == 1L, 1L] <- 1
  X[cbind(which(M == 1L), xi + 1L)] <- 1
  Q <- matrix(0, n, length(model$instrument_spec$levels) + 1L)
  Q[, 1L] <- 1; Q[cbind(seq_len(n), zi + 1L)] <- 1
  moments <- crossprod(Q, X) / n
  response <- colMeans(Q)
  covariance <- crossprod(Q) / n + model$tuning$weight_ridge * diag(ncol(Q))
  whitening <- t(chol(covariance))
  transformed_X <- forwardsolve(whitening, moments)
  transformed_y <- forwardsolve(whitening, response)
  penalty <- sqrt(model$tuning$lambda) * diag(ncol(X))[-1L, , drop = FALSE]
  augmented_X <- rbind(transformed_X, penalty)
  augmented_y <- c(transformed_y, rep(0, nrow(penalty)))
  coefficients <- qr.solve(augmented_X, augmented_y, tol = 1e-12)
  raw_package <- model$coefficients[1L] + model$coefficients[-1L]
  raw_reference <- coefficients[1L] + coefficients[-1L]
  deviation <- max(abs(raw_package - raw_reference))
  z <- copy(points[fold == j])
  z[, represented := cell %in% model$target_spec$levels]
  z[, history := sub('[|][^|]*[|][^|]*$', '', cell)]
  error_mass <- sum(z$probability * (z$estimate - z$true)^2)
  absent <- z[R2 == 1L & !represented]
  # Every unrepresented input has the same intercept in this fitted class.
  # Unequal true values on this set prove that it cannot contain the unique
  # population bridge on the complete finite support.
  distinct_absent <- uniqueN(round(absent$true, 12L))
  summary[[j]] <- data.table(fold = j, population_combinations = nrow(z),
    represented_combinations = sum(z$represented),
    unrepresented_combinations = sum(!z$represented),
    unrepresented_measured_probability = sum(z[represented == FALSE]$probability) / sum(z$probability),
    represented_fraction_of_squared_error = sum(z[represented == TRUE]$probability *
      (z[represented == TRUE]$estimate - z[represented == TRUE]$true)^2) / error_mass,
    true_values_on_unrepresented_R2_observed = distinct_absent,
    class_contains_unique_bridge_on_full_support = distinct_absent <= 1L)
  partition <- z[, .(combinations = .N,
    probability_conditional_measured = sum(probability) / sum(z$probability),
    rmse = sqrt(weighted.mean((estimate - true)^2, probability))), by = .(R2, represented)]
  partition[, fold := j]; partitions[[j]] <- partition
  for (history_value in unique(z[R2 == 1L]$history)) {
    ii <- rows[Hkeys[rows] == history_value]
    # Study the empirical conditional equation within each observed history.
    # Include all four reachable next-health categories, even if not sampled.
    aa <- action_index(d[ii, ], 2L, g)
    cc <- rep(NA_integer_, length(ii))
    measured <- d$R3[ii] == 1L
    cc[measured] <- 1L + d$C3_covariate[ii][measured] + 2L * d$Y3[ii][measured]
    counts <- tabulate(aa, nbins = nrow(g$A))
    J <- matrix(0, nrow(g$A), 4L)
    if (any(measured)) for (k in which(measured)) J[aa[k], cc[k]] <- J[aa[k], cc[k]] + 1
    active <- counts > 0
    if (any(active)) {
      operator <- J[active, , drop = FALSE] / counts[active]
      singular_values <- svd(operator, nu = 0L, nv = 0L)$d
      rank <- sum(singular_values > 1e-10 * max(c(singular_values, 1)))
      condition_number <- if (rank == 4L) max(singular_values) / min(singular_values) else Inf
    } else { rank <- 0L; condition_number <- Inf }
    zp <- z[history == history_value & R2 == 1L]
    histories[[length(histories) + 1L]] <- data.table(fold = j, history = history_value,
      training_n = length(ii), measured_n = sum(measured),
      sampled_treatment_combinations = sum(active),
      sampled_next_health_combinations = sum(colSums(J) > 0),
      empirical_equation_rank = rank, empirical_condition_number = condition_number,
      probability = sum(zp$probability),
      squared_error_mass = sum(zp$probability * (zp$estimate - zp$true)^2))
  }
  direct_checks[[j]] <- data.table(fold = j,
    max_raw_prediction_difference_from_independent_QR = deviation,
    objective_package = sum((augmented_X %*% model$coefficients - augmented_y)^2),
    objective_QR = sum((augmented_X %*% coefficients - augmented_y)^2))
}
summary <- rbindlist(summary); partitions <- rbindlist(partitions)
histories <- rbindlist(histories); direct_checks <- rbindlist(direct_checks)
fwrite(summary, file.path(out, 'actual-function-class-audit.csv'))
fwrite(partitions, file.path(out, 'prediction-error-by-support.csv'))
fwrite(histories, file.path(out, 'history-training-information.csv'))
fwrite(direct_checks, file.path(out, 'independent-sieve-calculation.csv'))
cat('Actual function class\n'); print(summary)
cat('Independent calculation\n'); print(direct_checks)
cat('Training information within R2=1 histories\n')
print(histories[, .(histories = .N, median_training_n = median(training_n),
  median_measured_n = median(measured_n),
  histories_rank_below_four = sum(empirical_equation_rank < 4L),
  median_full_rank_condition_number = median(empirical_condition_number[empirical_equation_rank == 4L]),
  fraction_of_error_mass_from_rank_below_four = sum(squared_error_mass[empirical_equation_rank < 4L]) / sum(squared_error_mass)), by = fold])
