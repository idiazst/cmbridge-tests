# Trace the points away from y=x using the saved million-person fit.
# Reconstruct its objective and optimality equations; no learner refitting.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
inputs <- readRDS(file.path(out, 'scoring-inputs.rds'))
fit <- readRDS(file.path(out, 'large-fit.rds'))
candidate <- fit$beta_fit$candidates[[1L]]
zkey <- cell_key(inputs$B); zlevels <- sort(unique(zkey))
zi <- match(zkey, zlevels)
xi <- match(inputs$beta_key[inputs$observed], candidate$target_levels)
stopifnot(!anyNA(xi))
counts <- tabulate(zi, length(zlevels))
joint_counts <- Matrix::sparseMatrix(i = zi[inputs$observed], j = xi,
  x = 1, dims = c(length(zlevels), length(candidate$target_levels)))
A <- Matrix::Diagonal(x = 1 / counts) %*% joint_counts
weights <- counts / sum(counts)
deviations <- candidate$coefficients[-1L]
beta <- candidate$coefficients[1L] + deviations
stopifnot(all(beta > candidate$tuning$lower & beta < candidate$tuning$upper))
moment_residual <- as.vector(A %*% beta) - 1
gradient <- as.vector(Matrix::crossprod(A, weights * moment_residual))
# glmnet rescales penalty factors to sum to the number of predictors.
effective_penalty <- candidate$tuning$penalty * (length(deviations) + 1L) / length(deviations)
zero <- abs(deviations) < 1e-9
kkt <- ifelse(zero, pmax(0, abs(gradient) - effective_penalty),
              abs(gradient + effective_penalty * sign(deviations)))
constant_gradient <- sum(weights * as.numeric(Matrix::rowSums(A)) * moment_residual)

points <- fread(file.path(out, 'bridge-large-errors.csv'))
index <- match(points$cell, candidate$target_levels)
stopifnot(!anyNA(index), max(abs(points$estimate - beta[index])) < 1e-10)
points[, `:=`(deviation_coefficient = deviations[index], zero_deviation = zero[index],
  objective_gradient = gradient[index], effective_penalty = effective_penalty,
  gradient_to_penalty_ratio = abs(gradient[index]) / effective_penalty,
  kkt_violation = kkt[index])]

saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
p <- saved$population
g <- make_mechanism('discrete_dose', dose_max = 3L)
pt <- study_task(p, g)$task
H <- pt$vars$history('A', 2L)
pop_observed <- p$R3 == 1L
V <- cbind(pt$natural[, H, drop = FALSE], p[, c('C3_covariate', 'Y3'), drop = FALSE])
key <- cell_key(V[pop_observed, , drop = FALSE])
first <- match(points$cell, key)
decoded <- p[which(pop_observed)[first], c('C1_baseline','A1_policy','A1_other','C2_covariate','Y2','C3_covariate','Y3')]
for (name in names(decoded)) set(points, j = name, value = decoded[[name]])

empirical_H <- cell_key(inputs$B[, H, drop = FALSE])
zhistory <- empirical_H[match(zlevels, zkey)]
target_history <- cell_key(V[pop_observed, H, drop = FALSE])[match(candidate$target_levels, key)]
operators <- list()
for (history in unique(empirical_H)) {
  ri <- which(zhistory == history); ci <- which(target_history == history)
  block <- as.matrix(A[ri, ci, drop = FALSE])
  spectrum <- svd(block)$d
  operators[[length(operators) + 1L]] <- data.table(history = history,
    training_people = sum(counts[ri]), measured = sum(joint_counts[ri, ]),
    minimum_treatment_group_n = min(counts[ri]), maximum_treatment_group_n = max(counts[ri]),
    operator_rank = sum(spectrum > max(spectrum) * 1e-10),
    operator_condition_number = max(spectrum) / min(spectrum))
}
operators <- rbindlist(operators)
points[, history := target_history[index]]
points <- merge(points, operators, by = 'history', all.x = TRUE, sort = FALSE)
fwrite(points[order(-abs(error))], file.path(out, 'large-bridge-point-diagnosis.csv'))
fwrite(operators, file.path(out, 'large-empirical-bridge-operators.csv'))

summary <- data.table(n = fit$n, training_n = fit$training_n,
  fitted_common_constant = candidate$coefficients[1L], selected_penalty = candidate$tuning$penalty,
  effective_penalty = effective_penalty, zero_deviation_coefficients = sum(zero),
  maximum_coefficient_kkt_violation = max(kkt), constant_gradient = constant_gradient,
  points_with_absolute_error_over_half = sum(abs(points$error) > .5),
  their_measured_population_probability = sum(points$probability[abs(points$error) > .5]),
  their_minimum_measured_training_n = min(points$training_measured[abs(points$error) > .5]),
  their_maximum_measured_training_n = max(points$training_measured[abs(points$error) > .5]),
  their_zero_deviation_coefficients = sum(points$zero_deviation[abs(points$error) > .5]),
  minimum_empirical_operator_rank = min(operators$operator_rank),
  maximum_empirical_operator_condition_number = max(operators$operator_condition_number))
fwrite(summary, file.path(out, 'large-bridge-point-summary.csv'))
print(summary, digits = 10)
print(points[order(-abs(error))][1:8, .(true,estimate,error,training_measured,
  probability,zero_deviation,gradient_to_penalty_ratio,minimum_treatment_group_n,
  operator_rank,operator_condition_number)], digits = 8)
cat('This is an algebraic check of the saved fitted objective, not a zero-penalty or known-function fit.\n')
