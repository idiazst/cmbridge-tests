# Same million-person binary sample and outer splits as the sieve bridge check.
# Only the saturated sieve adjoint is fitted. Its loading uses estimated,
# training-out-of-fold treatment ratios from the original SuperLearner library.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
library(ggplot2)
root <- Sys.getenv('SIM_OUTPUT')
input <- file.path(root, 'diagnostics/binary-dose-n1000000-seed5103006-sieve/large-sieve-fits.rds')
reference <- readRDS(input)
stopifnot(reference$n == 1000000L, reference$seed == 5103006L)
out <- file.path(root, 'diagnostics/binary-dose-n1000000-seed5103006-adjoint')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(out, 'large-adjoint-fits.rds')
if (file.exists(checkpoint)) stop('Completed checkpoint already exists; preserve it.')
g <- make_mechanism('binary_longitudinal', visits = 'missing_both')
require_valid_mechanism(g)
args <- study_arguments(g)
candidate <- cm_library(bridge = FALSE)[['sieve_md']]
stopifnot(candidate$method == 'sieve_md', candidate$control$target_basis == 'cell',
  candidate$control$instrument_basis == 'cell', is.null(candidate$control$lower),
  is.null(candidate$control$upper))
started <- proc.time()[['elapsed']]
cat(format(Sys.time()), 'Recreating the original million-person sample.\n')
set.seed(reference$seed)
d <- draw_data(reference$n, g)
stopifnot(identical(.Random.seed, reference$after_draw_seed))
prepared <- study_task(d, g, folds = reference$folds, learner_groups = 3L,
  balance_measurement = TRUE)
task <- prepared$task
stopifnot(identical(task$folds, reference$folds))
rm(prepared)
p <- reference$population
stopifnot(identical(p, enumerate_data(g)))
pt <- study_task(p, g)$task
H <- task$vars$history('A', 2L); A <- task$vars$A[[2L]]
B <- as.matrix(task$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(task$natural[, H, drop = FALSE], d[, args$health[[2L]], drop = FALSE]))
M <- d$R3
population_B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
population_V <- as.matrix(cbind(pt$natural[, H, drop = FALSE], p[, args$health[[2L]], drop = FALSE]))
bkey <- cell_key(population_B)
observed <- p$R3 == 1L
vkey <- cell_key(population_V[observed, , drop = FALSE])
hkey <- cell_key(pt$natural[, H, drop = FALSE])
# Verify uniqueness of the adjoint equation on all reachable histories.
uniqueness <- rbindlist(lapply(unique(hkey), function(history_value) {
  ii <- which(hkey == history_value & observed)
  z <- data.table(B = bkey[ii], V = vkey[match(ii, which(observed))],
    probability = p$probability[ii])[, .(probability = sum(probability)), by = .(B, V)]
  bl <- sort(unique(z$B)); vl <- sort(unique(z$V))
  operator <- matrix(0, length(vl), length(bl))
  operator[cbind(match(z$V, vl), match(z$B, bl))] <- z$probability
  operator <- operator / rowSums(operator)
  values <- svd(operator, nu = 0L, nv = 0L)$d
  rank <- sum(values > max(values) * 1e-10)
  data.table(history = history_value, R2 = p$R2[ii[1L]],
    reachable_treatment_values = length(bl), equation_rank = rank,
    unique_solution = rank == length(bl))
}))
stopifnot(all(uniqueness$unique_solution))
fwrite(uniqueness, file.path(out, 'population-uniqueness.csv'))
writeLines(c(capture.output(sessionInfo()),
  paste('bridge checkpoint SHA256', digest::digest(file = input, algo = 'sha256')),
  'n=1000000; seed=5103006; all original bridge outer splits reused.',
  'Saturated joint-category sieve adjoint; original positive ridge penalties.',
  'Ratio library: SL.saturated_l1_cv and SL.mean; positive penalties selected within training samples.',
  'Adjoint responses use treatment ratios fitted out of sample inside each outer training sample.'),
  file.path(out, 'session.txt'))
slim_sl <- function(fit) {
  fit[c('Y', 'Z', 'SL.predict', 'library.predict', 'validRows', 'call')] <- NULL
  fit
}
models <- ratio_models <- points <- equation_points <- metrics <- ratios_summary <- selections <- loading_records <- list()
for (j in seq_along(task$folds)) {
  rows <- task$folds[[j]]$training_set
  labels <- task$learner_folds[[j]]
  nested <- lmtp:::make_bridge_nested_folds(task$id[rows], d[rows, args$measurement, drop = FALSE], labels)
  control <- study_control(3L)
  control$.nested_split_function <- nested
  training_ratios <- matrix(NA_real_, length(rows), task$tau)
  inner_models <- list()
  population_loading <- rep(0, nrow(p))
  for (label in sort(unique(labels))) {
    local_train <- which(labels != label); local_valid <- which(labels == label)
    train <- rows[local_train]; valid <- rows[local_valid]
    stopifnot(!any(task$id[train] %in% task$id[valid]),
      !any(train %in% task$folds[[j]]$validation_set), !any(valid %in% task$folds[[j]]$validation_set))
    sub <- task$clone(deep = FALSE)
    sub$folds <- list(list(training_set = train, validation_set = valid))
    sub$learner_folds <- list(nested(task$id[train]))
    cat(format(Sys.time()), 'Outer fit', j, 'ratio response fold', label,
      ': training', length(train), 'evaluation', length(valid), '\n')
    ratio_started <- proc.time()[['elapsed']]
    fitted <- lmtp:::estimate_r(sub, 1L, trt_library, FALSE, control, function(...) NULL)
    stopifnot(!anyNA(fitted$ratios), all(is.finite(fitted$ratios)),
      !any(vapply(fitted$fits, function(f) any(f$errorsInCVLibrary) || any(f$errorsInLibrary), FALSE)))
    training_ratios[local_valid, ] <- fitted$ratios
    fitted$fits <- lapply(fitted$fits, slim_sl)
    inner_models[[as.character(label)]] <- fitted$fits
    population_ratios <- predict_ratios(fitted$fits, pt, g)
    population_loading[observed] <- population_loading[observed] + length(local_valid) / length(rows) *
      apply(population_ratios[observed, , drop = FALSE], 1L, prod) * p$Y3[observed]
    ww <- as.data.table(fit_candidate_weights(fitted$fits, 'ratio', j))
    ww[, response_fold := label]
    ss <- as.data.table(fit_tuning(fitted$fits, 'ratio', j))
    ss[, response_fold := label]
    selections[[length(selections) + 1L]] <- list(weights = ww, tuning = ss)
    truth_ratios <- population_truth_functions(p, g)$ratios
    ratios_summary[[length(ratios_summary) + 1L]] <- data.table(fold = j, response_fold = label,
      training_n = length(train), evaluation_n = length(valid),
      ratio1_rmse = sqrt(sum(p$probability * (population_ratios[, 1L] - truth_ratios[, 1L])^2)),
      ratio2_rmse = sqrt(sum(p$probability * (population_ratios[, 2L] - truth_ratios[, 2L])^2)),
      seconds = proc.time()[['elapsed']] - ratio_started)
    rm(fitted, sub)
    gc()
  }
  stopifnot(all(is.finite(training_ratios)))
  loading <- apply(training_ratios, 1L, prod) * ifelse(is.na(d$Y3[rows]), 0, d$Y3[rows])
  cat(format(Sys.time()), 'Fitting saturated sieve adjoint', j, 'on', sum(M[rows]), 'measured observations.\n')
  model <- cmbridge::fit_adjoint(B[rows, , drop = FALSE], V[rows, , drop = FALSE],
    M[rows], loading, method = candidate$method, control = candidate$control)
  stopifnot(model$tuning$lambda == 1e-8, model$tuning$weight_ridge == 1e-8)
  estimate <- predict(model, population_B)
  truth <- population_truth_functions(p, g)
  z <- data.table(cell = bkey, R2 = p$R2, true = truth$lambda[, 2L], estimate = estimate,
    probability = p$probability)[, {
      stopifnot(diff(range(true)) < 1e-10, diff(range(estimate)) < 1e-10, uniqueN(R2) == 1L)
      .(R2 = R2[1L], true = true[1L], estimate = estimate[1L], probability = sum(probability))
    }, by = cell]
  z[, `:=`(fold = j, represented = cell %in% model$target_spec$levels)]
  stopifnot(all(is.finite(z$estimate)))
  omega <- apply(truth$ratios, 1L, prod)
  target <- conditional_mean(omega[observed] * p$Y3[observed], p$probability[observed], vkey)
  left <- conditional_mean(estimate[observed], p$probability[observed], vkey)
  estimated_target <- conditional_mean(population_loading[observed], p$probability[observed], vkey)
  q <- data.table(cell = vkey, R2 = p$R2[observed], true = unname(target[vkey]),
    estimate = unname(left[vkey]), estimated_ratio_target = unname(estimated_target[vkey]),
    probability = p$probability[observed])[,
    .(R2 = R2[1L], true = true[1L], estimate = estimate[1L],
      estimated_ratio_target = estimated_ratio_target[1L], probability = sum(probability)), by = cell]
  q[, fold := j]
  for (r in c(1L, 0L)) {
    zz <- z[R2 == r]; qq <- q[R2 == r]
    metrics[[length(metrics) + 1L]] <- data.table(fold = j, R2 = r,
      training_n = length(rows), measured_training_n = sum(M[rows]),
      target_combinations = nrow(zz), represented_combinations = sum(zz$represented),
      adjoint_rmse = sqrt(weighted.mean((zz$estimate - zz$true)^2, zz$probability)),
      adjoint_signed_error = weighted.mean(zz$estimate - zz$true, zz$probability),
      adjoint_max_absolute_error = max(abs(zz$estimate - zz$true)),
      equation_rmse = sqrt(weighted.mean((qq$estimate - qq$true)^2, qq$probability)),
      equation_max_absolute_error = max(abs(qq$estimate - qq$true)),
      estimated_ratio_target_rmse = sqrt(weighted.mean((qq$estimated_ratio_target - qq$true)^2, qq$probability)))
  }
  beta <- reference$models[[j]]$predict_fun(population_V[observed, , drop = FALSE])
  bridge_residual <- rep(-1, nrow(p)); bridge_residual[observed] <- beta - 1
  bias <- sum(p$probability[observed] * omega[observed] * p$Y3[observed] * beta) - exact_truth(g)[2L] -
    sum(p$probability * estimate * bridge_residual)
  product <- -sum(p$probability[observed] * (estimate[observed] - truth$lambda[observed, 2L]) *
    (beta - truth$beta[observed, 2L]))
  stopifnot(abs(bias - product) < 1e-8)
  loading_records[[j]] <- data.table(fold = j, bridge_remainder_using_true_treatment_ratios = bias,
    product_identity_error = bias - product)
  keep <- c('method', 'type', 'coefficients', 'target_spec', 'instrument_spec',
    'target_design_cols', 'moment_loss', 'tuning', 'predict_fun', 'n', 'target_dim', 'instrument_dim')
  models[[j]] <- model[keep]
  class(models[[j]]) <- 'cmbridge_fit'
  ratio_models[[j]] <- inner_models
  points[[j]] <- z; equation_points[[j]] <- q
  fwrite(rbindlist(metrics), file.path(out, 'fit-diagnostics.csv'))
  fwrite(rbindlist(ratios_summary), file.path(out, 'ratio-fit-diagnostics.csv'))
  saveRDS(list(model = models[[j]], ratio_models = inner_models, fold = j,
    points = z, equation_points = q, loading_record = loading_records[[j]]),
    file.path(out, sprintf('completed-fold-%d.rds', j)))
  cat(format(Sys.time()), 'Completed outer adjoint fit', j, '\n')
  rm(model, training_ratios, inner_models)
  gc()
}
points <- rbindlist(points); equation_points <- rbindlist(equation_points)
weights <- rbindlist(lapply(selections, `[[`, 'weights'))
tuning <- rbindlist(lapply(selections, `[[`, 'tuning'))
fwrite(points, file.path(out, 'adjoint-predictions.csv'))
fwrite(equation_points, file.path(out, 'adjoint-equation-predictions.csv'))
fwrite(weights, file.path(out, 'ratio-learner-weights.csv'))
fwrite(tuning, file.path(out, 'ratio-penalty-selection.csv'))
fwrite(rbindlist(loading_records), file.path(out, 'bridge-remainder-diagnostic.csv'))
plot_points <- function(values, equation = FALSE) {
  values <- copy(values)
  values[, `:=`(fit = factor(paste('Training fit', fold), levels = paste('Training fit', 1:3)),
    measurement = factor(ifelse(R2 == 1L, 'R₂ = 1', 'R₂ = 0'), levels = c('R₂ = 1', 'R₂ = 0')))]
  limits <- range(values$true, values$estimate)
  limits <- limits + c(-1, 1) * max(.05, diff(limits) * .04)
  ggplot(values, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
    geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
    scale_size_area(max_size = 5, guide = 'none') +
    coord_equal(xlim = limits, ylim = limits, expand = FALSE) + facet_grid(measurement ~ fit) +
    labs(title = if (equation) 'Binary treatment: saturated sieve_md adjoint equation' else
      'Binary treatment: saturated sieve_md adjoint, n = 1,000,000',
      subtitle = 'Seed 5103006. Outcome at time 3. Same outer splits; estimated out-of-sample treatment ratios.',
      x = if (equation) 'True right-hand side of the adjoint equation' else expression('True '*lambda[2](H[2], A[2])),
      y = if (equation) 'Conditional mean of the fitted adjoint' else expression('Estimated '*lambda[2](H[2], A[2])),
      caption = paste('Blue line: y = x. Point area: population probability',
        if (equation) 'among measured final outcomes.' else 'of the history/treatment combination.',
        '\nThe binary adjoint is unique on the reachable support for both R₂ values.')) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
      plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 10),
      plot.margin = margin(12, 16, 12, 12))
}
pdf_device <- if (capabilities('aqua')) {
  function(filename, ...) grDevices::quartz(type = 'pdf', file = filename, ...)
} else grDevices::cairo_pdf
for (equation in c(FALSE, TRUE)) {
  figure <- plot_points(if (equation) equation_points else points, equation)
  name <- if (equation) 'adjoint-equation-true-versus-estimated' else 'adjoint-true-versus-estimated'
  ggsave(file.path(out, paste0(name, '.png')), figure, width = 11.5, height = 8.5, dpi = 190, bg = 'white')
  ggsave(file.path(out, paste0(name, '.pdf')), figure, width = 11.5, height = 8.5, device = pdf_device, bg = 'white')
}
saveRDS(list(models = models, ratio_models = ratio_models, folds = reference$folds,
  learner_folds = task$learner_folds, n = reference$n, seed = reference$seed,
  candidate = candidate, population = p, uniqueness = uniqueness,
  seconds = proc.time()[['elapsed']] - started), checkpoint)
cat('Adjoint recovery:\n'); print(rbindlist(metrics))
cat('Bridge remainder diagnostic:\n'); print(rbindlist(loading_records))
cat('Finished adjoint-only check in', proc.time()[['elapsed']] - started, 'seconds.\n')
