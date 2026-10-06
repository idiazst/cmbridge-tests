# Compare the original positive penalty candidates on the saved training data.
# No generating function is supplied to any fit; truth is evaluation only.
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

scenarios <- data.table(scale = c(.001,.01,.1,.01), tolerance = c(1e-8,1e-8,1e-8,1e-12),
  scenario = c('Original lower positive candidate','Original CV-selected candidate',
               'Original higher positive candidate','Selected candidate, tighter solver tolerance'))
results <- predictions <- models <- list()
for (i in seq_len(nrow(scenarios))) {
  control <- candidate$tuning
  control$penalty <- scenarios$scale[i] / sqrt(large$training_n)
  control$tolerance <- scenarios$tolerance[i]
  started <- proc.time()[3L]
  fitted <- cmbridge::fit_bridge(inputs$B, V, as.numeric(inputs$observed),
    method = 'saturated_l1', control = control)
  estimate <- predict(fitted, target)
  if (i == 2L) stopifnot(max(abs(estimate - points$estimate)) < 1e-8)
  tab <- copy(points)
  set(tab, j = 'estimate', value = estimate)
  tab[, `:=`(scenario = scenarios$scenario[i],
             scale = scenarios$scale[i], tolerance = scenarios$tolerance[i])]
  predictions[[i]] <- tab
  results[[i]] <- data.table(scenario = scenarios$scenario[i], scale = scenarios$scale[i],
    tolerance = scenarios$tolerance[i], penalty = control$penalty,
    original_cv_loss = large$tuning[kind == 'beta' & scale == scenarios$scale[i]]$loss,
    original_cv_selected = large$tuning[kind == 'beta' & scale == scenarios$scale[i]]$selected,
    beta_rmse = sqrt(weighted.mean((estimate - points$true)^2, points$probability)),
    maximum_absolute_error = max(abs(estimate - points$true)),
    points_with_absolute_error_over_half = sum(abs(estimate - points$true) > .5),
    their_measured_population_probability = sum(points$probability[abs(estimate - points$true) > .5]) / sum(points$probability),
    zero_deviation_coefficients = sum(abs(fitted$coefficients[-1L]) < 1e-9),
    fitted_constant = fitted$coefficients[1L],
    maximum_difference_from_original_selected = max(abs(estimate - points$estimate)),
    training_moment_loss = fitted$moment_loss, seconds = proc.time()[3L] - started)
  models[[i]] <- fitted[c('coefficients','target_levels','instrument_levels','tuning','moment_loss')]
}
fwrite(rbindlist(results), file.path(out, 'positive-penalty-shrinkage-comparison.csv'))
fwrite(rbindlist(predictions), file.path(out, 'positive-penalty-bridge-predictions.csv'))
saveRDS(list(scenarios = scenarios, models = models), file.path(out, 'positive-penalty-bridge-fits.rds'))
print(rbindlist(results), digits = 9)
cat('Only original positive penalty candidates were fitted on the same saved training sample.\n')
cat('The CV-selected fit was reproduced; no true functions or probabilities were used for fitting.\n')
