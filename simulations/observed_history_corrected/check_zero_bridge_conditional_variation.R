# Conditional sampling diagnostic for the history with the worst zero-L1 error.
# DGP probabilities generate categorical counts. The unpenalized calculation
# uses only those counts; true bridges are evaluation targets, never inputs.
library(data.table)
library(ggplot2)
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
tab <- fread(file.path(out, 'zero-penalty-worst-history-equations.csv'))
points <- fread(file.path(out, 'zero-penalty-point-diagnosis.csv'))
tab[, treatment := paste(A2_policy, A2_other, sep = '|')]
rows <- unique(tab$treatment); columns <- unique(tab$cell)
probability <- observed_counts <- matrix(0, length(rows), length(columns))
n <- integer(length(rows)); truth <- actual <- numeric(length(columns))
for (i in seq_len(nrow(tab))) {
  r <- match(tab$treatment[i], rows); c <- match(tab$cell[i], columns)
  probability[r, c] <- tab$population_joint_probability[i]
  observed_counts[r, c] <- tab$measured_health_n[i]
  n[r] <- tab$training_treatment_history_n[i]
  truth[c] <- tab$true_bridge[i]; actual[c] <- tab$unpenalized_bridge[i]
}
stopifnot(all(rowSums(probability) < 1), all(probability > 0),
          max(abs(probability %*% truth - 1)) < 1e-12)
weighted_A <- observed_counts / n * sqrt(n / sum(n))
actual_check <- qr.coef(qr(weighted_A, tol = 1e-12), sqrt(n / sum(n)))
stopifnot(max(abs(actual - actual_check)) < 1e-10)
selected_cell <- points$cell[1L]
target <- match(selected_cell, columns)
target_truth <- truth[target]; target_actual <- actual[target]
replications <- 10000L; seed <- 5204003L
set.seed(seed)
draws <- lapply(seq_along(n), function(row)
  rmultinom(replications, n[row], c(probability[row, ], 1 - sum(probability[row, ]))))
estimates <- matrix(NA_real_, replications, length(columns))
ranks <- integer(replications); condition <- numeric(replications)
for (b in seq_len(replications)) {
  joint <- do.call(rbind, lapply(draws, function(draw) draw[seq_along(columns), b]))
  A <- joint / n * sqrt(n / sum(n))
  factorization <- qr(A, tol = 1e-10)
  ranks[b] <- factorization$rank
  spectrum <- svd(A, nu = 0L, nv = 0L)$d
  condition[b] <- max(spectrum) / min(spectrum)
  if (ranks[b] == length(columns))
    estimates[b, ] <- qr.coef(factorization, sqrt(n / sum(n)))
}
valid <- ranks == length(columns)
records <- data.table(replication = seq_len(replications), full_rank = valid,
  empirical_weighted_condition_number = condition,
  target_unbounded_estimate = estimates[, target],
  target_bounded_estimate = pmax(1, pmin(6, estimates[, target])))
metrics <- rbindlist(lapply(c('No prediction bounds', 'Original bounds 1 to 6'), function(name) {
  beta <- if (name == 'No prediction bounds') records$target_unbounded_estimate[valid] else
    records$target_bounded_estimate[valid]
  data.table(scenario = name, seed = seed, replications = replications,
    full_rank_replications = sum(valid), rank_deficient_replications = sum(!valid),
    training_history_n = sum(n), minimum_treatment_group_n = min(n),
    maximum_treatment_group_n = max(n), true = target_truth,
    actual_estimate = target_actual, actual_error = target_actual - target_truth,
    sampling_mean = mean(beta), sampling_sd = sd(beta),
    sampling_rmse = sqrt(mean((beta - target_truth)^2)),
    quantile_025 = unname(quantile(beta, .025)),
    quantile_50 = unname(quantile(beta, .5)),
    quantile_975 = unname(quantile(beta, .975)),
    probability_absolute_error_at_least_observed = mean(abs(beta - target_truth) >= abs(target_actual - target_truth)),
    probability_estimate_at_least_observed = mean(beta >= target_actual))
}))
fwrite(metrics, file.path(out, 'zero-penalty-worst-history-sampling-summary.csv'))
fwrite(records, file.path(out, 'zero-penalty-worst-history-sampling-draws.csv'))
configurations <- data.table(cell = columns, true = truth,
  actual_unpenalized_estimate = actual,
  measured_training_n = colSums(observed_counts),
  expected_measured_training_n_conditional_on_history_counts = colSums(probability * n),
  sampling_mean = colMeans(estimates[valid, , drop = FALSE]),
  sampling_sd = apply(estimates[valid, , drop = FALSE], 2L, sd))
fwrite(configurations, file.path(out, 'zero-penalty-worst-history-sampling-health.csv'))
print(metrics, digits = 9)
print(configurations, digits = 9)
cat('Rank-deficient draws, if any, are reported separately and are not silently removed.\n')
cat('This is conditional sampling at the observed eight group sizes, not the full repeated SDR/TMLE study.\n')
cat('The target history was selected for its observed error: tail frequencies are descriptive, not a formal p-value.\n')
plot <- ggplot(records[full_rank == TRUE], aes(target_unbounded_estimate)) +
  geom_histogram(binwidth = .1, fill = '#718c82', color = 'white', linewidth = .15) +
  geom_vline(xintercept = target_truth, color = '#3265a8', linewidth = .8) +
  geom_vline(xintercept = target_actual, color = '#b85832', linewidth = .8) +
  labs(title = 'Sampling variation within the history with the largest unpenalized error',
    subtitle = sprintf('10,000 conditional samples; same treatment-group sizes, totaling %s people.', sum(n)),
    x = 'Unpenalized bridge estimate', y = 'Number of conditional samples',
    caption = sprintf('Blue: true bridge %.3f. Orange: actual estimate %.3f. History selected after examining its error.',
      target_truth, target_actual)) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(12, 12, 12, 12))
for (extension in c('png', 'pdf')) ggsave(file.path(out, 'figures', paste0('zero-penalty-conditional-sampling.', extension)),
  plot, width = 10, height = 5.5, dpi = 200, bg = 'white')
