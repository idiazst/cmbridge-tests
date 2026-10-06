# Algebraic diagnostics of saved fits, not additional estimator simulations.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/dose-ensemble-large')
new <- readRDS(file.path(out, 'large-fit.rds'))
old <- readRDS('results/observed-history-corrected/dose-cmbridge-debug/large-fit.rds')
inputs <- readRDS('results/observed-history-corrected/dose-cmbridge-debug/scoring-inputs.rds')
g <- make_mechanism('discrete_dose', dose_max = 3L)
p <- new$population; pt <- study_task(p, g)$task; truth <- population_truth_functions(p, g)
M <- p$R3 == 1L; H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(pt$natural[, H, drop = FALSE], p[, c('C3_covariate', 'Y3')]))
stopifnot(identical(colnames(B), colnames(inputs$B)), nrow(inputs$B) == new$training_n)
counts <- data.table(cell = inputs$beta_key[inputs$observed])[, .(training_measured = .N), by = cell]
history <- data.table(history = cell_key(inputs$B[, H, drop = FALSE]))[, .(history_training = .N), by = history]
pop_key <- cell_key(V[M, , drop = FALSE])
decoded <- as.data.table(p[M, c('C1_baseline', 'A1_policy', 'A1_other', 'C2_covariate', 'Y2', 'C3_covariate', 'Y3')])
decoded[, `:=`(cell = pop_key, history = cell_key(V[M, H, drop = FALSE]))]
decoded <- unique(decoded, by = 'cell')
support <- merge(merge(decoded, counts, by = 'cell', all.x = TRUE), history, by = 'history', all.x = TRUE)
support[is.na(training_measured), training_measured := 0L]
stopifnot(nrow(support) == 256L, !anyNA(support$history_training))
points <- rbindlist(lapply(c('bridge', 'adjoint-equation'), function(kind) {
  x <- fread(file.path(out, paste0(kind, '-predictions.csv')))
  x <- merge(x, support, by = 'cell', all.x = TRUE)
  x[, `:=`(function_name = kind, error = estimate - true)]
  x
}), use.names = TRUE)
fwrite(points, file.path(out, 'predictions-with-training-counts.csv'))
worst <- points[, .SD[order(-abs(error))][1:5], by = .(function_name, candidate)]
fwrite(worst, file.path(out, 'largest-point-errors.csv'))

features <- function(model, x) {
  if (model$method == 'saturated_l1') {
    index <- match(cell_key(x), model$target_levels)
    seen <- which(!is.na(index)); n <- nrow(x)
    return(as.matrix(Matrix::sparseMatrix(i = c(seq_len(n), seen),
      j = c(rep(1L, n), index[seen] + 1L), x = 1, dims = c(n, length(model$target_levels) + 1L))))
  }
  if (model$method == 'pmmr') return(cmbridge:::.nystrom_apply(
    cmbridge:::.scale_apply(x, model$target_scale), model$target_nystrom))
  cmbridge:::.eval_basis_spec(x, model$target_spec)
}
project <- function(A, target, weights) {
  weighted <- A * sqrt(weights); target <- target * sqrt(weights)
  sv <- svd(weighted)
  keep <- sv$d > max(sv$d) * 1e-10
  fit <- as.numeric(sv$u[, keep, drop = FALSE] %*% crossprod(sv$u[, keep, drop = FALSE], target))
  list(rank = sum(keep), rmse = sqrt(sum((fit - target)^2) / sum(weights)),
       max_residual = max(abs((fit - target) / sqrt(weights))))
}
conditional_matrix <- function(features, keys, weights) {
  levels <- sort(unique(keys)); index <- match(keys, levels)
  probs <- as.numeric(Matrix::sparseMatrix(i = index, j = rep(1L, length(index)), x = weights,
    dims = c(length(levels), 1L)))
  sums <- Matrix::sparseMatrix(i = index, j = seq_along(index), x = weights,
    dims = c(length(levels), length(index))) %*% features
  list(matrix = as.matrix(sums) / probs, probability = probs, levels = levels)
}
bkeys <- cell_key(B); blevels <- sort(unique(bkeys)); bi <- match(bkeys, blevels)
bp <- as.numeric(Matrix::sparseMatrix(i = bi, j = rep(1L, nrow(p)), x = p$probability,
                                     dims = c(length(blevels), 1L)))
omega <- apply(truth$ratios, 1L, prod)
phi <- conditional_mean(omega[M] * p$Y3[M], p$probability[M], pop_key)
classes <- list()
for (i in seq_along(new$candidate)) {
  name <- names(new$candidate)[i]
  beta_model <- new$beta_fit$candidates[[i]]
  lambda_model <- new$adjoint_fit$candidates[[i]]
  F <- features(beta_model, V[M, , drop = FALSE])
  span <- project(F, truth$beta[M, 2L], p$probability[M])
  beta_A <- as.matrix(Matrix::sparseMatrix(i = bi[M], j = seq_len(sum(M)), x = p$probability[M],
    dims = c(length(blevels), sum(M))) %*% F) / bp
  beta_equation <- project(beta_A, rep(1, length(bp)), bp)
  L <- features(lambda_model, B[M, , drop = FALSE])
  lambda_A <- conditional_matrix(L, pop_key, p$probability[M])
  lambda_equation <- project(lambda_A$matrix, unname(phi[lambda_A$levels]), lambda_A$probability)
  classes[[name]] <- data.table(candidate = name, beta_basis_columns = ncol(F), beta_basis_rank = span$rank,
    closest_true_beta_rmse_in_class = span$rmse, smallest_bridge_equation_rmse_in_class = beta_equation$rmse,
    adjoint_basis_columns = ncol(L), adjoint_equation_rank = lambda_equation$rank,
    smallest_adjoint_equation_rmse_in_class = lambda_equation$rmse,
    contains_true_bridge = span$rmse < 1e-8, contains_valid_adjoint_solution = lambda_equation$rmse < 1e-8)
}
classes <- rbindlist(classes)
fwrite(classes, file.path(out, 'function-class-audit.csv'))
verification <- data.table(
  maximum_saturated_beta_difference_from_previous_fit = max(abs(new$candidate$saturated_cv$beta - old$beta[, 2L])),
  maximum_saturated_lambda_difference_from_previous_fit = max(abs(new$candidate$saturated_cv$lambda - old$lambda[, 2L])),
  minimum_measured_target_cell_n = min(support$training_measured),
  measured_target_cells_under_10 = sum(support$training_measured < 10L),
  measured_target_cells_under_100 = sum(support$training_measured < 100L))
stopifnot(verification$maximum_saturated_beta_difference_from_previous_fit < 1e-10,
          verification$maximum_saturated_lambda_difference_from_previous_fit < 1e-10,
          classes[candidate %in% c('sieve_md', 'saturated_cv'), all(contains_true_bridge & contains_valid_adjoint_solution)])
fwrite(verification, file.path(out, 'verification.csv'))
print(classes, digits = 9); print(verification, digits = 9)
print(worst[candidate == 'ensemble', .(function_name, true, estimate, error, training_measured,
                                      history_training, probability)], digits = 9)
cat('Class containment is an algebraic audit only; no true values were used to fit the saved estimators.\n')
