# User-requested zero-L1 diagnostic on the SAME large training sample.
# No generating function or probability enters fitting.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
inputs <- readRDS(file.path(out, 'scoring-inputs.rds'))
large <- readRDS(file.path(out, 'large-fit.rds'))
candidate <- large$beta_fit$candidates[[1L]]
points <- fread(file.path(out, 'bridge-predictions.csv'))[source == 'large check']
keys <- unique(inputs$beta_key)
parsed <- do.call(rbind, strsplit(keys, '|', fixed = TRUE))
values <- matrix(NA_real_, nrow(parsed), ncol(parsed))
known <- parsed != 'NA'
values[known] <- as.numeric(parsed[known])
V <- values[match(inputs$beta_key, keys), , drop = FALSE]
target <- values[match(points$cell, keys), , drop = FALSE]
stopifnot(all(is.finite(V[inputs$observed, ])), all(is.finite(target)))

control <- candidate$tuning
control$penalty <- 0
control$tolerance <- 1e-12
started <- proc.time()[3L]
warnings <- character()
fitted <- withCallingHandlers(cmbridge::fit_bridge(inputs$B, V,
  as.numeric(inputs$observed), method = 'saturated_l1', control = control),
  warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
bounded <- predict(fitted, target)
index <- match(points$cell, fitted$target_levels)
stopifnot(!anyNA(index))
raw <- fitted$coefficients[1L] + fitted$coefficients[1L + index]

# Independent unpenalized weighted least-squares solution of the same
# empirical conditional-moment criterion, using only training counts.
zkey <- cell_key(inputs$B); zlevels <- sort(unique(zkey))
zi <- match(zkey, zlevels)
xi <- match(inputs$beta_key[inputs$observed], fitted$target_levels)
counts <- tabulate(zi, length(zlevels))
joint <- Matrix::sparseMatrix(i = zi[inputs$observed], j = xi, x = 1,
  dims = c(length(zlevels), length(fitted$target_levels)))
A <- as.matrix(Matrix::Diagonal(x = 1 / counts) %*% joint)
weighted_A <- A * sqrt(counts / sum(counts))
weighted_y <- sqrt(counts / sum(counts))
direct <- qr.coef(qr(weighted_A, tol = 1e-12), weighted_y)
stopifnot(all(is.finite(direct)))
direct_raw <- direct[index]
direct_bounded <- pmax(control$lower, pmin(control$upper, direct_raw))

saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
p <- saved$population
g <- make_mechanism('discrete_dose', dose_max = 3L)
pt <- study_task(p, g)$task
truth <- population_truth_functions(p, g)
H <- pt$vars$history('A', 2L)
M <- p$R3 == 1L
pop_key <- cell_key(cbind(pt$natural[M, H, drop = FALSE],
                         p[M, c('C3_covariate','Y3'), drop = FALSE]))
pop_index <- match(pop_key, points$cell)
stopifnot(!anyNA(pop_index))
omega <- apply(truth$ratios, 1L, prod)
theta <- exact_truth(g)[2L]
scenarios <- list('CV-selected positive penalty, original bounds' = points$estimate,
                  'Zero L1 penalty, original bounds' = bounded,
                  'Zero L1 penalty, no prediction bounds' = raw,
                  'Direct zero-penalty solution, original bounds' = direct_bounded,
                  'Direct zero-penalty solution, no prediction bounds' = direct_raw)
summary <- tables <- list()
for (name in names(scenarios)) {
  estimate <- scenarios[[name]]
  tab <- copy(points); set(tab, j = 'estimate', value = estimate)
  tab[, scenario := name]
  tables[[length(tables) + 1L]] <- tab
  residual <- rep(-1, nrow(p)); residual[M] <- estimate[pop_index] - 1
  implied_error <- sum(p$probability[M] * omega[M] * p$Y3[M] * estimate[pop_index]) - theta
  correction <- -sum(p$probability * large$lambda[,2L] * residual)
  original_far <- abs(points$estimate - points$true) > .5
  summary[[length(summary) + 1L]] <- data.table(scenario = name,
    penalty = if (grepl('CV-selected', name)) candidate$tuning$penalty else 0,
    beta_rmse = sqrt(weighted.mean((estimate - points$true)^2, points$probability)),
    mean_prediction_error = weighted.mean(estimate - points$true, points$probability),
    maximum_absolute_error = max(abs(estimate - points$true)),
    minimum_prediction = min(estimate), maximum_prediction = max(estimate),
    points_with_absolute_error_over_half = sum(abs(estimate - points$true) > .5),
    improved_among_original_16_far_points = sum(original_far & abs(estimate - points$true) < abs(points$estimate - points$true)),
    beta_representation_error = implied_error,
    bridge_remainder_with_unchanged_adjoint = implied_error + correction)
}
diagnostic <- data.table(training_n = large$training_n,
  empirical_operator_rank = qr(weighted_A)$rank, target_coefficients = ncol(weighted_A),
  maximum_package_vs_direct_raw_prediction_difference = max(abs(raw - direct_raw)),
  predictions_clipped_to_lower_bound = sum(raw < control$lower),
  predictions_clipped_to_upper_bound = sum(raw > control$upper),
  zero_fit_seconds = proc.time()[3L] - started,
  warnings = paste(unique(warnings), collapse = ' | '))
fwrite(rbindlist(summary), file.path(out, 'zero-penalty-bridge-comparison.csv'))
fwrite(rbindlist(tables), file.path(out, 'zero-penalty-bridge-predictions.csv'))
fwrite(diagnostic, file.path(out, 'zero-penalty-solver-check.csv'))
saveRDS(list(coefficients = fitted$coefficients, target_levels = fitted$target_levels,
  tuning = fitted$tuning, direct_coefficients = direct, diagnostic = diagnostic),
  file.path(out, 'zero-penalty-bridge-fit.rds'))
print(rbindlist(summary), digits = 9)
print(diagnostic, digits = 9)
cat('True functions are used only for evaluation; both zero-penalty solutions use empirical training counts.\n')
cat('The adjoint is unchanged. This is a bridge component check, not a final SDR/TMLE or coverage study.\n')
