# Exact representability checks; no fitting to known functions or new simulation.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
g <- make_mechanism('discrete_dose', dose_max = 3L)
p <- saved$population
pt <- study_task(p, g)$task
nu <- population_truth_functions(p, g)
H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
M <- p$R3 == 1L
B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(pt$natural[, H, drop = FALSE], p[, c('C3_covariate', 'Y3'), drop = FALSE]))
bkey <- cell_key(B); vkey <- cell_key(V[M, , drop = FALSE])
bs <- sort(unique(bkey)); vs <- sort(unique(vkey))
joint <- as.matrix(Matrix::sparseMatrix(i = match(vkey, vs), j = match(bkey[M], bs),
  x = p$probability[M], dims = c(length(vs), length(bs))))
vprob <- rowSums(joint)
Tstar <- joint / vprob
phi <- unname(conditional_mean(apply(nu$ratios[M, , drop = FALSE], 1L, prod) * p$Y3[M],
                              p$probability[M], vkey)[vs])
beta_true <- vapply(vs, function(cell) unique(nu$beta[M, 2L][vkey == cell]), 0)
lambda_true <- vapply(bs, function(cell) unique(nu$lambda[, 2L][bkey == cell]), 0)
stopifnot(max(abs(Tstar %*% lambda_true - phi)) < 1e-10,
          all(beta_true >= 1 & beta_true <= 6))

results <- list()
append_beta <- function(levels, n, fold, source) {
  absent <- !(vs %in% levels)
  error <- rep(0, length(vs))
  if (any(absent)) {
    constant <- weighted.mean(beta_true[absent], vprob[absent])
    error[absent] <- constant - beta_true[absent]
  }
  data.table(n = n, fold = fold, source = source, function_name = 'beta',
    population_targets = length(vs), represented_targets = sum(!absent),
    missing_target_probability = sum(vprob[absent]) / sum(vprob),
    exact_member_or_solution = max(abs(error)) < 1e-10,
    maximum_representation_or_equation_error = max(abs(error)),
    minimum_probability_weighted_error = sqrt(weighted.mean(error^2, vprob)))
}
append_adjoint <- function(levels, n, fold, source) {
  seen <- bs %in% levels
  # The actual fitted dictionary permits arbitrary deviations on seen
  # treatment/history cells and only one common constant on all unseen cells.
  D <- cbind(1, diag(length(bs))[, seen, drop = FALSE])
  E <- Tstar %*% D
  decomposition <- svd(E)
  keep <- decomposition$d > max(decomposition$d) * 1e-10
  coefficients <- decomposition$v[, keep, drop = FALSE] %*%
    (as.vector(crossprod(decomposition$u[, keep, drop = FALSE], phi)) / decomposition$d[keep])
  solution <- as.vector(D %*% coefficients)
  residual <- as.vector(Tstar %*% solution - phi)
  fwrite(data.table(cell = bs, value = solution), file.path(out,
    sprintf('class-adjoint-solution-n%d-fold%d.csv', n, fold)))
  data.table(n = n, fold = fold, source = source, function_name = 'adjoint',
    population_targets = length(bs), represented_targets = sum(seen),
    missing_target_probability = sum(p$probability[!bkey %in% levels]),
    exact_member_or_solution = max(abs(residual)) < 1e-8,
    maximum_representation_or_equation_error = max(abs(residual)),
    minimum_probability_weighted_error = NA_real_)
}
results[[1L]] <- append_beta(vs, 0L, 0L, 'Complete population dictionary')
results[[2L]] <- append_adjoint(bs, 0L, 0L, 'Complete population dictionary')

set.seed(saved$design$seed)
d <- draw_data(saved$design$n, g)
prepared <- study_task(d, g, folds = saved$folds, learner_folds = saved$learner_folds,
                       learner_groups = 3L)
dt <- prepared$task
for (j in seq_along(saved$folds)) {
  rows <- saved$folds[[j]]$training_set
  observed <- d$R3[rows] == 1L
  target_beta <- cell_key(as.matrix(cbind(dt$natural[rows, H, drop = FALSE],
                           d[rows, c('C3_covariate', 'Y3'), drop = FALSE])))
  target_adjoint <- cell_key(as.matrix(dt$natural[rows, c(H, A), drop = FALSE]))
  results[[length(results) + 1L]] <- append_beta(unique(target_beta[observed]), nrow(d), j, 'Realized training dictionary')
  results[[length(results) + 1L]] <- append_adjoint(unique(target_adjoint[observed]), nrow(d), j, 'Realized training dictionary')
}
large <- readRDS(file.path(out, 'large-fit.rds'))
results[[length(results) + 1L]] <- append_beta(large$beta_fit$candidates[[1L]]$target_levels,
                                            large$n, large$fold, 'Realized large training dictionary')
results[[length(results) + 1L]] <- append_adjoint(large$adjoint_fit$candidates[[1L]]$target_levels,
                                               large$n, large$fold, 'Realized large training dictionary')
answer <- rbindlist(results)
fwrite(answer, file.path(out, 'function-class-membership.csv'))
print(answer, digits = 8)
cat('Positive penalties affect the fitted coefficients, not this algebraic span check.\n')
cat('The code constructs its dictionaries only from observed training target cells.\n')
cat('For the adjoint, membership means at least one valid solution, not equality to a chosen solution.\n')
