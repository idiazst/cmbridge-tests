source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
library(ggplot2)
out <- 'results/observed-history-corrected/dose-ensemble-r002'
new <- readRDS(file.path(out, 'discrete_dose-n4000-r002.rds'))
old <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
stopifnot(nrow(new$errors) == 0L, identical(old$population, new$population), identical(old$folds, new$folds),
          identical(old$learner_folds, new$learner_folds))
p <- new$population; g <- make_mechanism('discrete_dose', dose_max = 3L)
pt <- study_task(p, g)$task; truth <- population_truth_functions(p, g)
M <- p$R3 == 1L; H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
B <- cell_key(pt$natural[, c(H, A), drop = FALSE])
V <- cell_key(cbind(pt$natural[M, H, drop = FALSE], p[M, c('C3_covariate', 'Y3')]))
omega <- apply(truth$ratios, 1L, prod); theta <- exact_truth(g)[2L]
phi <- conditional_mean(omega[M] * p$Y3[M], p$probability[M], V)
metrics <- beta_points <- adjoint_points <- list()
for (name in c('Previous saturated L1 only', 'Expanded library')) {
  result <- if (name == 'Expanded library') new else old
  for (fold in seq_along(result$folds)) {
    nu <- result$fitted_functions[[paste('sdr', 'beta1_lambda1_ratio1_m1', fold, sep = '/')]]
    residual <- rep(-1, nrow(p)); residual[M] <- nu$beta[M, 2L] - 1
    bridge_equation <- conditional_mean(residual, p$probability, B)
    adjoint_equation <- conditional_mean(nu$lambda[M, 2L], p$probability[M], V)
    beta_error <- sum(p$probability[M] * omega[M] * p$Y3[M] * nu$beta[M, 2L]) - theta
    correction <- -sum(p$probability * nu$lambda[, 2L] * residual)
    metrics[[length(metrics) + 1L]] <- data.table(library = name, fold = fold,
      beta_rmse = sqrt(weighted.mean((nu$beta[M, 2L] - truth$beta[M, 2L])^2, p$probability[M])),
      beta_representation_error = beta_error, bridge_correction = correction,
      bridge_remainder = beta_error + correction,
      bridge_equation_rmse = sqrt(sum(p$probability * unname(bridge_equation[B])^2)),
      adjoint_equation_rmse = sqrt(weighted.mean((unname(adjoint_equation[V]) - unname(phi[V]))^2,
                                                 p$probability[M])))
    bp <- data.table(cell = V, true = truth$beta[M, 2L], estimate = nu$beta[M, 2L],
      probability = p$probability[M])[, .(true = unique(true), estimate = unique(estimate),
      probability = sum(probability)), by = cell]
    bp[, `:=`(library = name, fold = fold)]
    beta_points[[length(beta_points) + 1L]] <- bp
    ap <- bp[, .(cell, probability, library, fold)]
    ap[, `:=`(true = unname(phi[cell]), estimate = unname(adjoint_equation[cell]))]
    adjoint_points[[length(adjoint_points) + 1L]] <- ap
  }
}
metrics <- rbindlist(metrics); beta_points <- rbindlist(beta_points)
adjoint_points <- rbindlist(adjoint_points)
fwrite(metrics, file.path(out, 'component-comparison.csv'))
fwrite(beta_points, file.path(out, 'bridge-predictions.csv'))
fwrite(adjoint_points, file.path(out, 'adjoint-equation-predictions.csv'))
dir.create(file.path(out, 'figures'), showWarnings = FALSE)
plot_comparison <- function(points, title, file, xlab, ylab, bounded = FALSE) {
  points <- copy(points)
  points[, library := factor(library, levels = c('Previous saturated L1 only', 'Expanded library'))]
  limit <- range(points$true, points$estimate)
  limit <- limit + c(-1, 1) * diff(limit) * .04
  figure <- ggplot(points, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .65) +
    geom_point(aes(size = probability), color = '#254a40', alpha = .6) +
    scale_size_area(max_size = 4, guide = 'none') +
    coord_equal(xlim = limit, ylim = limit, expand = FALSE) + facet_grid(library ~ fold) +
    labs(title = title, subtitle = 'Original n = 4,000 sample, seed 5203002; original estimator and learner splits.',
      x = xlab, y = ylab,
      caption = 'Columns: outer training samples 1, 2, 3. Blue: y = x. Point area: probability among measured outcomes.') +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.spacing = unit(1.1, 'lines'),
      plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 9))
  for (extension in c('png', 'pdf')) ggsave(file.path(out, 'figures', paste0(file, '.', extension)),
    figure, width = 11, height = 8, dpi = 180, bg = 'white')
}
plot_comparison(beta_points, 'Bridge at time 2: previous and expanded libraries', 'bridge-comparison',
                'True bridge from the DGP', 'Estimated bridge', TRUE)
plot_comparison(adjoint_points, 'Adjoint equation at time 2: previous and expanded libraries', 'adjoint-equation-comparison',
                'True right-hand side of adjoint equation', 'Left-hand side evaluated over the DGP')
weights <- new$ensemble_weights
cm <- weights[kind %in% c('beta', 'adjoint')]
cm_table <- dcast(cm, candidate ~ kind + fold, value.var = 'weight')
fwrite(cm_table, file.path(out, 'cmbridge-weight-comparison.csv'))
regression_table <- dcast(weights[kind == 'regression' & fold == 2L],
  estimator + horizon_or_depth ~ method, value.var = 'weight')
fwrite(regression_table, file.path(out, 'regression-weights-sample2.csv'))
stopifnot(!any(weights[kind == 'ratio', method] == 'SL.earth'),
  all(weights[, abs(sum(weight) - 1) < 1e-8, by = .(kind, fold, horizon_or_depth, estimator)]$V1))
cm[, function_name := ifelse(kind == 'beta', 'Bridge', 'Adjoint')]
figure <- ggplot(cm, aes(candidate, weight, fill = candidate)) + geom_col() +
  facet_grid(function_name ~ fold) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'cmbridge candidate weights in the n = 4,000 rerun',
    subtitle = 'Weights selected from shared learner-validation samples.', x = NULL, y = 'Ensemble weight') +
  theme_minimal(base_size = 11) + theme(legend.position = 'none',
    axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = 'bold'))
for (extension in c('png', 'pdf')) ggsave(file.path(out, 'figures', paste0('cmbridge-weights.', extension)),
  figure, width = 11, height = 6, dpi = 180, bg = 'white')

comparison <- fread(file.path(out, 'old-versus-expanded-library.csv'))[horizon == 2L]
selected_old <- metrics[library == 'Previous saturated L1 only' & fold == 2L]
selected_new <- metrics[library == 'Expanded library' & fold == 2L]
lines <- c('# Numerical dose example: original n = 4,000 sample with the expanded learner libraries', '',
  'This rerun uses seed 5203002 and the original estimator and learner assignments. Both sets of assignments match the saved run exactly. The two treatment times, dose values 0–3, data-generating mechanism, and true mean are unchanged. The true mean at outcome time 3 is 0.871138516760.', '',
  '## Learners and package changes', '',
  'For both bridge and adjoint, cmbridge combines `sieve_md` with saturated joint-category target and conditioning bases, `landweber` with its default additive spline bases, `pmmr` with its default Gaussian-kernel approximation, and `saturated_l1` with cross-validated positive penalties. Sieve and Landweber were timed using the same joint-category bases on training sample 2: median fitting times were 0.085 and 0.111 seconds, respectively. The faster sieve method receives the saturated basis in this run.', '',
  'The outcome regressions use `SL.saturated_l1_cv`, MARS (`SL.earth`), and `SL.mean`. Treatment-ratio classification uses `SL.saturated_l1_cv` and `SL.mean`. MARS is absent from all cmbridge and treatment-ratio calls. All candidate weights and validation losses are retained in [ensemble_weights.csv](ensemble_weights.csv). The detailed settings and splitting construction are in the [saved simulation description](source/reports/simulation-plan.md).', '',
  'The saturated basis includes every interaction among the full input variables. Its training dictionary has independent deviations only for combinations observed in that training sample; unseen target combinations receive the common constant. This preserves the distinction between the full saturated function class and finite training samples with empty combinations.', '',
  'cmbridge 0.3.0.9003 adds the joint-category basis and applies the specified bridge bounds to every candidate. lmtp remains 1.6.0.9004; earth is 5.3.6. The first expanded-library attempt failed because default B-splines could produce NaN predictions when a training-constant predictor varied in validation. The corrected basis retains a constant contribution for that unidentified effect. The failed checkpoint, log, and source are preserved in [attempt-1](attempt-1/run.log). Package tests covering interactions, bounds, and this spline case pass; see [cmbridge-tests.log](cmbridge-tests.log).', '',
  '## Final estimates at outcome time 3', '',
  '| Estimator | Previous estimate | Expanded-library estimate | Standard error | 95% interval |',
  '|:--|--:|--:|--:|:--|')
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i]
  lines <- c(lines, sprintf('| %s | %.6f | %.6f | %.6f | [%.6f, %.6f] |',
    toupper(row$estimator), row$old_estimate, row$estimate, row$se, row$lower, row$upper))
}
lines <- c(lines, '', 'Both 95% intervals at outcome time 3 still exclude the true value. These are results from one dataset; they do not estimate confidence-interval coverage or repeated-sample bias.', '',
  '## Bridge and adjoint diagnostics', '',
  'The following comparisons evaluate saved fits over the exact DGP distribution; true functions are not supplied to any learner. Training sample 2 was the most problematic in the previous run. The bridge contribution to the remainder isolates this component using true treatment ratios for evaluation, holding out the sequential-regression contribution.', '',
  '| Quantity, training sample 2 | Previous saturated L1 only | Expanded library |',
  '|:--|--:|--:|')
for (item in list(c('Bridge prediction RMSE','beta_rmse'),
  c('Bridge conditional-equation RMSE','bridge_equation_rmse'),
  c('Adjoint conditional-equation RMSE','adjoint_equation_rmse'),
  c('Bridge contribution to remainder','bridge_remainder'))) {
  lines <- c(lines, sprintf('| %s | %.6f | %.6f |', item[1L], selected_old[[item[2L]]], selected_new[[item[2L]]]))
}
lines <- c(lines, '', 'The bridge itself improves little: in sample 2, saturated L1 receives 98.56% weight and saturated sieve 1.44%. The improvement in the bridge correction comes mainly from the adjoint, whose equation error falls substantially. Both the remaining bridge error and the final interval misses are retained.', '',
  '![Bridge estimate versus truth](figures/bridge-comparison.png)', '',
  'The adjoint plot compares the two sides of its conditional equation. It does not demand equality to one arbitrarily selected adjoint solution.', '',
  '![Adjoint conditional equation](figures/adjoint-equation-comparison.png)', '',
  '[All sample-specific component comparisons](component-comparison.csv) retain each of the three training samples.', '',
  '## Candidate weights', '',
  '| Candidate | Bridge sample 1 | Bridge sample 2 | Bridge sample 3 | Adjoint sample 1 | Adjoint sample 2 | Adjoint sample 3 |',
  '|:--|--:|--:|--:|--:|--:|--:|')
for (i in seq_len(nrow(cm_table))) {
  row <- cm_table[i]
  lines <- c(lines, sprintf('| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |', row$candidate,
    row$beta_1, row$beta_2, row$beta_3, row$adjoint_1, row$adjoint_2, row$adjoint_3))
}
lines <- c(lines, '', '![cmbridge candidate weights](figures/cmbridge-weights.png)', '',
  'Outcome-regression weights in training sample 2 are:', '',
  '| Estimator | Pooled diagonal | Saturated L1 | MARS | Mean |', '|:--|--:|--:|--:|--:|')
for (i in seq_len(nrow(regression_table))) {
  row <- regression_table[i]
  lines <- c(lines, sprintf('| %s | %s | %.4f | %.4f | %.4f |', toupper(row$estimator),
    row$horizon_or_depth, row$SL.saturated_l1, row$SL.earth, row$SL.mean))
}
lines <- c(lines, '', 'Diagonal 1 pools m₂,₁ and m₃,₂; diagonal 2 fits m₃,₁. [The complete weight table](ensemble_weights.csv) includes all regression and treatment-ratio fits across all three training samples. A metadata assignment initially copied the candidate method into the SDR/TMLE label of regression-weight rows. [The repair](weight-label-repair.log) verifies every corrected label against independent saved tuning records; no fits, weights, losses, or estimates change. The original executed checkpoint is preserved.', '',
  'The run retains a warning about moving repeated interior spline knots inside the boundary. All candidates nevertheless completed with finite predictions; there are no candidate failures in the completed run.', '',
  sprintf('Both estimators completed in %.1f minutes. No failed candidate is silently assigned zero weight. The full repeated-sample study has not been launched with this new configuration.', new$seconds / 60), '',
  '[Numerical results](results.csv), [execution log](run.log), [reproduction instructions](REPRODUCE.md), and [frozen source checksums](SOURCE-SHA256SUMS) are included.')
writeLines(lines, file.path(out, 'REPORT.md'))
print(metrics, digits = 8)
print(cm_table, digits = 7)
print(regression_table, digits = 7)
cat('Warnings recorded by the rerun:', paste(new$job_warnings, collapse = ' | '), '\n')
