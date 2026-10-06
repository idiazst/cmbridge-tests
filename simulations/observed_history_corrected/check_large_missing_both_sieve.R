# Large-sample check of one existing sieve candidate, with no other learners.
# The frozen DGP, full-history predictors, regularization, and bounds are used.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
library(ggplot2)
root <- Sys.getenv('SIM_OUTPUT')
original <- readRDS(file.path(root, 'jobs/binary_longitudinal-n4000-r006.rds'))
n <- 1000000L
seed <- original$design$seed
stopifnot(seed == 5103006L, original$design$visits == 'missing_both')
out <- file.path(root, sprintf('diagnostics/binary-dose-n%d-seed%d-sieve', n, seed))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(out, 'large-sieve-fits.rds')
if (file.exists(checkpoint)) stop('Completed checkpoint already exists; preserve it.')
g <- make_mechanism('binary_longitudinal', visits = original$design$visits)
require_valid_mechanism(g)
p <- original$population
stopifnot(identical(p, enumerate_data(g)))
pt <- study_task(p, g)$task
args <- study_arguments(g)
H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
# Build only the matrices needed for the bridge. This duplicates the study's
# missing-value encoding, not restrictions on how the true function depends
# on the history. Verify it against the original lmtp task below.
bridge_predictors <- function(d) {
  e <- as.data.frame(d)
  for (t in seq_len(g$tau)) for (name in args$time_vary[[t]]) {
    key <- paste0('..bridge_', t, '_', name)
    e[[key]] <- ifelse(is.na(d[[name]]), 0, d[[name]])
    e[[paste0(key, '_missing')]] <- as.numeric(is.na(d[[name]]))
  }
  list(B = as.matrix(e[, c(H, A), drop = FALSE]),
    V = as.matrix(cbind(e[, H, drop = FALSE], d[, args$health[[2L]], drop = FALSE])))
}
reference <- readRDS(file.path(root, 'diagnostics/binary-dose-n4000-r006-sieve/recovered-sieve-fits.rds'))
small_task <- study_task(reference$data, g, folds = original$folds,
  learner_folds = original$learner_folds, learner_groups = original$design$learner_groups,
  balance_measurement = original$design$balance_measurement)$task
small <- bridge_predictors(reference$data)
stopifnot(identical(small$B, as.matrix(small_task$natural[, c(H, A), drop = FALSE])),
  identical(small$V, as.matrix(cbind(small_task$natural[, H, drop = FALSE],
    reference$data[, args$health[[2L]], drop = FALSE]))))
rm(reference, small_task, small)
gc()
candidate <- cm_library()[['sieve_md']]
stopifnot(candidate$method == 'sieve_md', candidate$control$target_basis == 'cell',
  candidate$control$instrument_basis == 'cell')
writeLines(c(capture.output(sessionInfo()),
  paste('n', n), paste('seed', seed), paste('RNG', paste(RNGkind(), collapse = ', ')),
  'Only saturated sieve_md for the outcome at time 3 is fitted.',
  'Full history and the original frozen library controls are retained.'), file.path(out, 'session.txt'))
cat(format(Sys.time()), 'Drawing', n, 'binary-treatment observations.\n')
set.seed(seed)
started <- proc.time()[['elapsed']]
d <- draw_data(n, g)
after_draw_seed <- .Random.seed
# Same frozen outer-split implementation as LmtpTask for distinct person IDs
# and continuous outcome_type; the bridge-only fit has no inner tuning step.
folds <- lmtp:::make_folds(data.frame(person = seq_len(n)), V = 3L, cluster_ids = seq_len(n))
pred <- bridge_predictors(d)
B <- pred$B; V <- pred$V
rm(pred)
Hkeys <- cell_key(B[, H, drop = FALSE])
M <- d$R3
health_category <- rep(NA_integer_, n)
health_category[M == 1L] <- 1L + d$C3_covariate[M == 1L] + 2L * d$Y3[M == 1L]
treatment_category <- action_index(d, 2L, g)
observed <- p$R3 == 1L
population_pred <- bridge_predictors(p)
newV <- population_pred$V[observed, , drop = FALSE]
population_B <- cell_key(population_pred$B)
key <- cell_key(newV)
truth <- population_truth_functions(p, g)$beta[observed, 2L]
stopifnot(all(truth >= candidate$control$lower), all(truth <= candidate$control$upper))
models <- points <- metrics <- histories <- list()
for (j in seq_along(folds)) {
  train <- folds[[j]]$training_set
  cat(format(Sys.time()), 'Fitting sieve bridge', j, 'on', length(train), 'training observations.\n')
  fit_started <- proc.time()[['elapsed']]
  model <- cmbridge::fit_bridge(B[train, , drop = FALSE], V[train, , drop = FALSE],
    M[train], method = candidate$method, control = candidate$control)
  fit_seconds <- proc.time()[['elapsed']] - fit_started
  stopifnot(model$tuning$lambda == 1e-8, model$tuning$weight_ridge == 1e-8)
  estimate <- predict(model, newV)
  z <- data.table(cell = key, R2 = p$R2[observed], true = truth, estimate = estimate,
    probability = p$probability[observed])[, {
      stopifnot(diff(range(true)) < 1e-12, diff(range(estimate)) < 1e-12, uniqueN(R2) == 1L)
      .(R2 = R2[1L], true = true[1L], estimate = estimate[1L], probability = sum(probability))
    }, by = cell]
  z[, `:=`(fold = j, represented = cell %in% model$target_spec$levels)]
  stopifnot(all(is.finite(z$estimate)))
  residual <- rep(-1, nrow(p)); residual[observed] <- estimate - 1
  moments <- conditional_mean(residual, p$probability, population_B)
  for (r in c(1L, 0L)) {
    zz <- z[R2 == r]
    metrics[[length(metrics) + 1L]] <- data.table(fold = j, R2 = r,
      n = n, training_n = length(train), measured_training_n = sum(M[train]),
      population_combinations = nrow(zz), represented_combinations = sum(zz$represented),
      rmse = sqrt(weighted.mean((zz$estimate - zz$true)^2, zz$probability)),
      signed_error = weighted.mean(zz$estimate - zz$true, zz$probability),
      max_absolute_error = max(abs(zz$estimate - zz$true)),
      bridge_equation_max_error = max(abs(moments)),
      target_cells = length(model$target_spec$levels), instrument_cells = length(model$instrument_spec$levels),
      fit_seconds = fit_seconds)
  }
  for (history_value in unique(Hkeys[train][d$R2[train] == 1L])) {
    ii <- train[Hkeys[train] == history_value]
    counts <- tabulate(treatment_category[ii], nbins = nrow(g$A))
    J <- matrix(0, nrow(g$A), 4L)
    selected <- ii[M[ii] == 1L]
    tab <- data.table(action = treatment_category[selected], health = health_category[selected])[
      , .(count = .N), by = .(action, health)]
    J[cbind(tab$action, tab$health)] <- tab$count
    active <- counts > 0
    operator <- J[active, , drop = FALSE] / counts[active]
    values <- svd(operator, nu = 0L, nv = 0L)$d
    rank <- sum(values > 1e-10 * max(c(values, 1)))
    histories[[length(histories) + 1L]] <- data.table(fold = j, history = history_value,
      training_n = length(ii), measured_n = sum(M[ii]), empirical_equation_rank = rank,
      empirical_condition_number = if (rank == 4L) max(values) / min(values) else Inf)
  }
  keep <- c('method', 'type', 'coefficients', 'target_spec', 'instrument_spec',
    'target_design_cols', 'moment_loss', 'tuning', 'predict_fun', 'n', 'target_dim', 'instrument_dim')
  models[[j]] <- model[keep]
  points[[j]] <- z
  rm(model)
  gc()
}
points <- rbindlist(points); metrics <- rbindlist(metrics); histories <- rbindlist(histories)
points[, `:=`(fit_label = factor(paste('Training fit', fold), levels = paste('Training fit', 1:3)),
  measurement_label = factor(ifelse(R2 == 1L, 'R₂ = 1', 'R₂ = 0'), levels = c('R₂ = 1', 'R₂ = 0')))]
plot_sieve <- function(values, all_measurement_states = FALSE) {
  limits <- range(values$true, values$estimate)
  limits <- limits + c(-1, 1) * max(.025, diff(limits) * .04)
  figure <- ggplot(values, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
    geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
    scale_size_area(max_size = 5, guide = 'none') +
    coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    labs(title = 'Binary treatment: saturated sieve_md bridge, n = 1,000,000',
      subtitle = paste('Seed 5103006. Bridge for outcome at time 3.',
        if (!all_measurement_states) 'R₂ = 1: the bridge solution is unique.' else ''),
      x = expression('True '*beta[2](H[2], C[3])),
      y = expression('Estimated '*beta[2](H[2], C[3])),
      caption = paste('Blue line: y = x. Point area: population probability among measured final outcomes.',
        if (all_measurement_states) 'R₂ = 0: the reference is one valid solution; it is not unique.' else '', sep = '\n')) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
      plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 10),
      plot.margin = margin(12, 16, 12, 12))
  if (all_measurement_states) figure + facet_grid(measurement_label ~ fit_label)
  else figure + facet_wrap(~ fit_label, nrow = 1)
}
pdf_device <- if (capabilities('aqua')) {
  function(filename, ...) grDevices::quartz(type = 'pdf', file = filename, ...)
} else grDevices::cairo_pdf
for (all_states in c(FALSE, TRUE)) {
  figure <- plot_sieve(if (all_states) points else points[R2 == 1L], all_states)
  name <- if (all_states) 'sieve-true-versus-estimated-all' else 'sieve-true-versus-estimated'
  height <- if (all_states) 8.5 else 4.8
  ggsave(file.path(out, paste0(name, '.png')), figure, width = 11.5, height = height, dpi = 190, bg = 'white')
  ggsave(file.path(out, paste0(name, '.pdf')), figure, width = 11.5, height = height, device = pdf_device, bg = 'white')
}
fwrite(points, file.path(out, 'bridge-predictions.csv'))
fwrite(metrics, file.path(out, 'fit-diagnostics.csv'))
fwrite(histories, file.path(out, 'history-training-information.csv'))
saveRDS(list(models = models, folds = folds, n = n, seed = seed, after_draw_seed = after_draw_seed,
  candidate = candidate, history_columns = H, treatment_columns = A, population = p,
  seconds = proc.time()[['elapsed']] - started), checkpoint)
cat('Bridge recovery:\n'); print(metrics)
cat('Training information at observed time 2:\n')
print(histories[, .(histories = .N, median_training_n = median(training_n),
  median_measured_n = median(measured_n), minimum_measured_n = min(measured_n),
  histories_rank_below_four = sum(empirical_equation_rank < 4L)), by = fold])
cat('Finished bridge-only check in', proc.time()[['elapsed']] - started, 'seconds.\n')
