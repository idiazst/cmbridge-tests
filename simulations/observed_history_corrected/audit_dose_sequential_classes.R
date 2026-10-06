# Check the actual sequential predictor spaces and classifier link.
# Exact functions are used only to test representability, never for fitting.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
saved <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
g <- make_mechanism('discrete_dose', dose_max = 3L)
p <- saved$population
pt <- study_task(p, g)$task
nu <- population_truth_functions(p, g)
keys <- targets <- list()
for (s in 1:2) {
  H <- pt$vars$history('A', s); A <- pt$vars$A[[s]]
  keys[[s]] <- cell_key(pt$natural[, c(H, A), drop = FALSE])
  terminal <- numeric(nrow(p)); measured <- p[[paste0('R', s + 1L)]] == 1L
  terminal[measured] <- p[[paste0('Y', s + 1L)]][measured] * nu$beta[measured, s]
  targets[[s]] <- conditional_mean(terminal, p$probability, keys[[s]])
}
H2 <- pt$vars$history('A', 2L); A2 <- pt$vars$A[[2L]]
shifted2 <- pt$natural[, c(H2, A2), drop = FALSE]
shifted2[, A2] <- pt$shifted[, A2, drop = FALSE]
earlier <- conditional_mean(unname(targets[[2L]][cell_key(shifted2)]), p$probability, keys[[1L]])
earlier_shift <- pt$natural[, c(pt$vars$history('A', 1L), pt$vars$A[[1L]]), drop = FALSE]
earlier_shift[, pt$vars$A[[1L]]] <- pt$shifted[, pt$vars$A[[1L]], drop = FALSE]
stopifnot(abs(sum(p$probability * unname(earlier[cell_key(earlier_shift)])) - exact_truth(g)[2L]) < 1e-12)

long <- lmtp:::pivot(pt$natural, pt$vars)
true_values <- list(c(unname(targets[[1L]][keys[[1L]]]), unname(targets[[2L]][keys[[2L]]])) / 6,
                    unname(earlier[keys[[1L]]]) / 6)
weights <- list(rep(p$probability / 2, 2L), p$probability)
features <- function(data, depth) {
  rows <- data$time <= 3L - depth
  vars <- setdiff(names(data), c('..i..C_1', '..i..C_1_lag', '..i..wide_id', '..i..N', '..i..D_1', '..i..R'))
  if (depth == 2L) vars <- setdiff(vars, 'time')
  for (name in setdiff(vars, c('..i..lmtp_id', '..i..Y_1'))) {
    if (all(duplicated(data[rows, name])[-1L])) vars <- setdiff(vars, name)
  }
  setdiff(vars, c('..i..lmtp_id', '..i..Y_1'))
}
set.seed(saved$design$seed)
d <- draw_data(saved$design$n, g)
dt <- study_task(d, g, folds = saved$folds, learner_folds = saved$learner_folds, learner_groups = 3L)$task
records <- list()
for (j in c(0L, 1:3)) for (depth in 1:2) {
  train_long <- if (j == 0L) long else lmtp:::pivot(dt$natural[saved$folds[[j]]$training_set, ], dt$vars)
  cols <- features(train_long, depth)
  poprows <- long$time <= 3L - depth
  trainrows <- train_long$time <= 3L - depth
  popkey <- lmtp:::.saturated_keys(long[poprows, cols, drop = FALSE])
  seen <- unique(lmtp:::.saturated_keys(train_long[trainrows, cols, drop = FALSE]))
  values <- true_values[[depth]]
  cells <- data.table(cell = popkey, truth = values, probability = weights[[depth]])[,
    .(truth_range = diff(range(truth)), truth = weighted.mean(truth, probability),
      probability = sum(probability)), by = cell]
  stopifnot(max(cells$truth_range) < 1e-12)
  absent <- !(cells$cell %in% seen)
  error <- numeric(nrow(cells))
  if (any(absent)) error[absent] <- weighted.mean(cells$truth[absent], cells$probability[absent]) - cells$truth[absent]
  records[[length(records) + 1L]] <- data.table(estimators = 'SDR and TMLE', fold = j,
    scope = if (j == 0L) 'Complete population dictionary' else 'Realized training dictionary',
    recursion = if (depth == 1L) 'm_2,1 and m_3,2, pooled by time' else 'm_3,1',
    population_predictor_combinations = nrow(cells), represented_combinations = sum(!absent),
    missing_predictor_probability = sum(cells$probability[absent]),
    exact_class_member = max(abs(error)) < 1e-10,
    minimal_scaled_prediction_rmse = sqrt(weighted.mean(error^2, cells$probability)),
    scaled_true_min = min(values), scaled_true_max = max(values),
    tmle_bounds_include_truth = all(values > study_control()$.bound & values < 1 - study_control()$.bound))
}
fwrite(rbindlist(records), file.path(out, 'sequential-regression-class-membership.csv'))
print(rbindlist(records), digits = 8)

# Equal natural/shifted stacking makes the classifier probability g_d/(g+g_d).
classification <- list()
for (s in 1:2) {
  state <- state_index(p, s); action <- action_index(p, s, g)
  probability <- g$gd[cbind(p$C1_baseline + 1L, s, state, action)] /
    (g$g[cbind(p$C1_baseline + 1L, s, state, action)] + g$gd[cbind(p$C1_baseline + 1L, s, state, action)])
  target_columns <- c(pt$vars$history('A', s), pt$vars$A[[s]])
  popkey <- lmtp:::.saturated_keys(pt$natural[, target_columns, drop = FALSE])
  cells <- data.table(cell = popkey, probability = probability, weight = p$probability)[,
    .(truth_range = diff(range(probability)), truth = weighted.mean(probability, weight), weight = sum(weight)), by = cell]
  stopifnot(max(cells$truth_range) < 1e-12)
  for (j in c(0L, 1:3)) {
    train <- if (j == 0L) pt$natural else dt$natural[saved$folds[[j]]$training_set, ]
    shifted <- if (j == 0L) pt$shifted else dt$shifted[saved$folds[[j]]$training_set, ]
    stack <- lmtp:::stack_data(train, shifted, dt$vars$A, dt$vars$C, s)
    seen <- unique(lmtp:::.saturated_keys(stack[, target_columns, drop = FALSE]))
    absent <- !cells$cell %in% seen
    absent_range <- if (any(absent)) diff(range(cells$truth[absent])) else 0
    classification[[length(classification) + 1L]] <- data.table(time = s, fold = j,
      scope = if (j == 0L) 'Complete population dictionary' else 'Realized training dictionary',
      population_combinations = nrow(cells), represented_combinations = sum(!absent),
      missing_natural_probability = sum(cells$weight[absent]),
      zero_probability_combinations = sum(cells$truth == 0),
      true_probability_min = min(cells$truth), true_probability_max = max(cells$truth),
      exact_finite_logistic_representation = all(cells$truth > 0 & cells$truth < 1) && absent_range < 1e-10,
      member_of_logistic_class_closure = absent_range < 1e-10,
      ratio_upper_cap_inactive_at_truth = all(cells$truth < .999))
  }
}
fwrite(rbindlist(classification), file.path(out, 'treatment-classifier-class-membership.csv'))
print(rbindlist(classification), digits = 8)
cat('The Gaussian saturated candidate includes the true sequential functions on complete support.\n')
cat('The logistic saturated candidate represents zero classifier probabilities only as a limit.\n')
cat('Training-only dictionaries also impose a common prediction on unseen joint combinations.\n')
