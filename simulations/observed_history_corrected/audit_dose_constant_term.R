# Algebraic diagnostic of the fitted learner's constant term, not a new fit.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
x <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
g <- make_mechanism('discrete_dose', dose_max = 3L)
set.seed(x$design$seed)
d <- draw_data(x$design$n, g)
prepared <- study_task(d, g, folds = x$folds, learner_folds = x$learner_folds,
                       learner_groups = 3L)
rows <- x$folds[[2L]]$training_set
H <- prepared$task$vars$history('A', 2L)
A <- prepared$task$vars$A[[2L]]
key <- cell_key(prepared$task$natural[rows, c(H, A), drop = FALSE])
counts <- data.table(cell = key, M = d$R3[rows])[, .(n = .N, q = mean(M)), by = cell]
empirical_mean_q <- weighted.mean(counts$q, counts$n)
empirical_mean_q2 <- weighted.mean(counts$q^2, counts$n)
p <- x$population
pt <- study_task(p, g)$task
popkey <- cell_key(pt$natural[, c(H, A), drop = FALSE])
q <- conditional_mean(p$R3, p$probability, popkey)
population_constant <- sum(p$probability * unname(q[popkey])) /
                       sum(p$probability * unname(q[popkey])^2)
predictions <- fread(file.path(out, 'bridge-reference-predictions.csv'))
most_common <- as.numeric(names(sort(table(predictions$estimate), decreasing = TRUE))[1L])
result <- data.table(
  training_n = length(rows), measured_fraction = empirical_mean_q,
  empirical_mean_squared_conditional_measurement = empirical_mean_q2,
  empirical_constant_criterion_minimizer = empirical_mean_q / empirical_mean_q2,
  population_constant_criterion_minimizer = population_constant,
  most_common_fitted_bridge_value = most_common)
fwrite(result, file.path(out, 'constant-term-criterion-diagnostic.csv'))
print(result, digits = 10)
cat('This evaluates the constant term of the existing criterion. It does not fit another learner.\n')
