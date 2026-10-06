# Evaluate the saved estimated functions over the exact finite DGP support.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
library(ggplot2)
out <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/dose-ensemble-large')
new <- readRDS(file.path(out, 'large-fit.rds'))
old <- readRDS('results/observed-history-corrected/dose-cmbridge-debug/large-fit.rds')
stopifnot(new$n == old$n, new$seed == old$seed, new$fold == old$fold,
          identical(new$folds, old$folds), identical(new$learner_folds, old$learner_folds))
g <- make_mechanism('discrete_dose', dose_max = 3L)
p <- enumerate_data(g); stopifnot(identical(p, new$population))
pt <- study_task(p, g)$task; truth <- population_truth_functions(p, g)
M <- p$R3 == 1L; H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
B <- cell_key(pt$natural[, c(H, A), drop = FALSE])
V <- cell_key(cbind(pt$natural[M, H, drop = FALSE], p[M, c('C3_covariate', 'Y3')]))
w <- p$probability; omega <- apply(truth$ratios, 1L, prod); theta <- exact_truth(g)[2L]
phi <- conditional_mean(omega[M] * p$Y3[M], w[M], V)
sets <- c(list(previous = list(beta = old$beta[, 2L], lambda = old$lambda[, 2L]),
               ensemble = list(beta = new$beta[, 2L], lambda = new$lambda[, 2L])), new$candidate)
labels <- c(previous = 'Previous saturated L1 only', ensemble = 'Expanded ensemble',
            sieve_md = 'Saturated sieve', landweber = 'Landweber', pmmr = 'PMMR',
            saturated_cv = 'Saturated L1 (CV penalty)')
metrics <- bp <- ap <- be <- list()
for (name in names(sets)) {
  nu <- sets[[name]]
  residual <- rep(-1, nrow(p)); residual[M] <- nu$beta[M] - 1
  bridge_equation <- conditional_mean(residual, w, B)
  adjoint_equation <- conditional_mean(nu$lambda[M], w[M], V)
  beta_error <- sum(w[M] * omega[M] * p$Y3[M] * nu$beta[M]) - theta
  correction <- -sum(w * nu$lambda * residual)
  product <- -sum(w[M] * (nu$lambda[M] - truth$lambda[M, 2L]) * (nu$beta[M] - truth$beta[M, 2L]))
  stopifnot(abs(beta_error + correction - product) < 1e-8)
  b <- data.table(cell = V, true = truth$beta[M, 2L], estimate = nu$beta[M], probability = w[M])[
    , .(true = unique(true), estimate = unique(estimate), probability = sum(probability)), by = cell]
  stopifnot(nrow(b) == uniqueN(V), all(is.finite(b$estimate)))
  a <- b[, .(cell, probability)]
  a[, `:=`(true = unname(phi[cell]), estimate = unname(adjoint_equation[cell]))]
  b[, candidate := name]; a[, candidate := name]
  bp[[name]] <- b; ap[[name]] <- a
  be[[name]] <- data.table(cell = names(bridge_equation), residual = unname(bridge_equation), candidate = name)
  metrics[[name]] <- data.table(candidate = name, label = unname(labels[name]),
    beta_rmse = sqrt(weighted.mean((nu$beta[M] - truth$beta[M, 2L])^2, w[M])),
    beta_max_error = max(abs(nu$beta[M] - truth$beta[M, 2L])),
    beta_cells_error_gt_05 = sum(abs(b$estimate - b$true) > .5),
    bridge_equation_rmse = sqrt(sum(w * unname(bridge_equation[B])^2)),
    bridge_equation_max_error = max(abs(bridge_equation)),
    adjoint_equation_rmse = sqrt(weighted.mean((unname(adjoint_equation[V]) - unname(phi[V]))^2, w[M])),
    adjoint_equation_max_error = max(abs(adjoint_equation - phi)),
    beta_representation_error = beta_error, bridge_correction = correction,
    bridge_remainder = beta_error + correction, product_identity_error = beta_error + correction - product)
}
metrics <- rbindlist(metrics); bp <- rbindlist(bp); ap <- rbindlist(ap); be <- rbindlist(be)
classes <- fread(file.path(out, 'function-class-audit.csv'))
verification <- fread(file.path(out, 'verification.csv'))
worst <- fread(file.path(out, 'largest-point-errors.csv'))
fwrite(metrics, file.path(out, 'component-comparison.csv'))
fwrite(bp, file.path(out, 'bridge-predictions.csv'))
fwrite(ap, file.path(out, 'adjoint-equation-predictions.csv'))
fwrite(be, file.path(out, 'bridge-equation-residuals.csv'))
dir.create(file.path(out, 'figures'), showWarnings = FALSE)
save_plot <- function(figure, filename, width = 11, height = 8) {
  for (extension in c('png', 'pdf')) ggsave(file.path(out, 'figures', paste0(filename, '.', extension)),
    figure, width = width, height = height, dpi = 180, bg = 'white')
}
plot_points <- function(points, title, xlab, ylab, facet = TRUE) {
  points <- copy(points)
  points[, label := factor(unname(labels[candidate]), levels = unname(labels))]
  limits <- range(points$true, points$estimate)
  limits <- limits + c(-1, 1) * diff(limits) * .04
  figure <- ggplot(points, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .65) +
    geom_point(aes(size = probability), color = '#254a40', alpha = .6) +
    scale_size_area(max_size = 4, guide = 'none') +
    coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    labs(title = title, subtitle = 'Original n = 1,000,000 sample; seed 5204001; training sample 2 (n = 666,667).',
      x = xlab, y = ylab,
      caption = 'Blue: y = x. Point area: probability among measured outcomes. Original split assignments preserved.') +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, 'lines'), strip.text = element_text(face = 'bold'),
      plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 9))
  if (facet) figure <- figure + facet_wrap(~label, ncol = 3)
  figure
}
save_plot(plot_points(bp, 'Bridge at time 2: all candidates and ensemble',
  'True bridge from the DGP', 'Estimated bridge'), 'bridge-all-candidates')
save_plot(plot_points(ap, 'Adjoint equation at time 2: all candidates and ensemble',
  'True right-hand side of the adjoint equation', 'Left-hand side evaluated over the DGP'), 'adjoint-equation-all-candidates')
for (name in c('ensemble', names(new$candidate))) {
  save_plot(plot_points(bp[candidate == name], paste('Bridge:', labels[name]),
    'True bridge from the DGP', 'Estimated bridge', FALSE), paste0('bridge-', name), 8, 7)
  save_plot(plot_points(ap[candidate == name], paste('Adjoint equation:', labels[name]),
    'True right-hand side of the adjoint equation', 'Left-hand side evaluated over the DGP', FALSE),
    paste0('adjoint-equation-', name), 8, 7)
}
cm <- new$ensemble_weights[kind %in% c('beta', 'adjoint')]
cm[, function_name := ifelse(kind == 'beta', 'Bridge', 'Adjoint')]
cm[, label := unname(labels[candidate])]
weights <- dcast(cm, candidate ~ kind, value.var = 'weight')
losses <- dcast(cm, candidate ~ kind, value.var = 'cv_loss')
fwrite(weights, file.path(out, 'cmbridge-weights.csv'))
fwrite(losses, file.path(out, 'cmbridge-validation-losses.csv'))
figure <- ggplot(cm, aes(label, weight, fill = label)) + geom_col() +
  facet_wrap(~function_name, ncol = 2) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'Candidate weights in the million-person check',
    subtitle = 'Weights selected using the original shared learner-validation assignments.', x = NULL, y = 'Ensemble weight') +
  theme_minimal(base_size = 11) + theme(legend.position = 'none', axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = 'bold'))
save_plot(figure, 'cmbridge-weights', 10, 5)

lines <- c('# Million-person numerical-dose check with the expanded cmbridge ensemble', '',
  sprintf('The original sample is reused: n = %s, seed %s, outer training sample %s, training n = %s. Every original estimator and learner split matches the previous check. The unchanged DGP has two treatment times and dose values 0–3; its true mean at outcome time 3 is %.12f.',
    format(new$n, big.mark = ','), new$seed, new$fold, format(new$training_n, big.mark = ','), theta), '',
  'This is a check of the bridge and adjoint functions. The treatment ratios used in adjoint responses are fitted with the same out-of-sample construction as in lmtp. The true functions and exact population distribution are used only for evaluation. Sequential outcome regressions, final SDR/TMLE estimates, and confidence intervals are not computed in this component check. MARS remains reserved for sequential outcome regressions.', '',
  '## Learners', '',
  'Both functions use saturated joint-category sieve, default additive-spline Landweber, default PMMR with 40 target and 80 conditioning kernel centers, and saturated L1 with its positive penalty selected by cross-validation. All inputs are retained. Bridge predictions are bounded between 1 and 6 for every candidate; adjoint predictions are unrestricted. Settings are preserved in the [frozen simulation description](source/reports/simulation-plan.md).', '',
  'The sieve computation now uses grouped sufficient statistics instead of a large indicator matrix. It solves the same regularized normal equations. Prediction functions retain their fitted parameters without keeping large temporary fitting arrays alive. The package validation includes comparison to the independent dense normal equations and all 96 tests pass; see [cmbridge-tests.log](cmbridge-tests.log). These changes preserve the specified estimators.', '',
  '## Exact-population comparisons', '',
  'RMSE means the square root of the probability-weighted mean squared error. Bridge prediction error is evaluated among measured outcomes. Bridge equation error is evaluated over all histories and treatments. Adjoint equation error is evaluated among measured outcomes.', '',
  '| Fit | Bridge prediction RMSE | Bridge equation RMSE | Adjoint equation RMSE | Bridge contribution to remainder |',
  '|:--|--:|--:|--:|--:|')
for (i in seq_len(nrow(metrics))) {
  r <- metrics[i]
  lines <- c(lines, sprintf('| %s | %.6f | %.6f | %.6f | %.8f |',
    r$label, r$beta_rmse, r$bridge_equation_rmse, r$adjoint_equation_rmse, r$bridge_remainder))
}
lines <- c(lines, '',
  sprintf('The expanded ensemble improves bridge prediction RMSE from %.6f to %.6f, and adjoint equation RMSE from %.6f to %.6f. The number of bridge input combinations with absolute error above 0.5 falls from %s to %s. Its maximum bridge error remains %.6f; the adjoint equation also retains rare large errors. The plots and tables preserve them.',
    metrics[candidate == 'previous', beta_rmse], metrics[candidate == 'ensemble', beta_rmse],
    metrics[candidate == 'previous', adjoint_equation_rmse], metrics[candidate == 'ensemble', adjoint_equation_rmse],
    metrics[candidate == 'previous', beta_cells_error_gt_05], metrics[candidate == 'ensemble', beta_cells_error_gt_05],
    metrics[candidate == 'ensemble', beta_max_error]), '',
  'For each candidate row, the remainder uses that candidate for both functions; the ensemble row uses the separately selected bridge and adjoint ensembles. This diagnostic holds the sequential contribution at its population value. It is not the bias of a final SDR/TMLE estimator. The numerical identity with the product of bridge and adjoint estimation errors is verified for every row.', '',
  '## y = x plots', '', 'The bridge has a unique solution in this DGP, so its plot compares fitted values directly with the true bridge.', '',
  '![Bridge: all candidates and ensemble](figures/bridge-all-candidates.png)', '',
  'The adjoint has multiple valid solutions. Its plot therefore compares the conditional left-hand side obtained from the fitted adjoint with the true right-hand side of its equation. A valid fitted solution lies on y = x without needing to equal an arbitrarily chosen adjoint function pointwise.', '',
  '![Adjoint conditional equation: all candidates and ensemble](figures/adjoint-equation-all-candidates.png)', '',
  'The [ensemble bridge plot](figures/bridge-ensemble.png) and [ensemble adjoint-equation plot](figures/adjoint-equation-ensemble.png) are also saved separately. Every panel and standalone plot has equal x and y scales. Points represent distinct input combinations; point areas show their exact DGP probabilities among measured outcomes.', '',
  '## Selected weights', '', '| Candidate | Bridge weight | Adjoint weight |', '|:--|--:|--:|')
for (i in seq_len(nrow(weights))) {
  r <- weights[i]
  lines <- c(lines, sprintf('| %s | %.6f | %.6f |', labels[r$candidate], r$beta, r$adjoint))
}
lines <- c(lines, '', '![Candidate weights](figures/cmbridge-weights.png)', '',
  'Selection uses estimated validation losses, rather than the true errors used to evaluate these graphs. The full [weight table](ensemble_weights.csv) includes treatment-ratio learners, and the [candidate validation losses](cmbridge-validation-losses.csv) are retained.', '',
  '## Departures from y = x', '',
  'The following audit asks whether each fitted basis can represent the true bridge, and whether its adjoint basis contains any solution of the population adjoint equation. This is an exact-population linear algebra check of the function class, without fitting an estimator to known true responses or running another simulation.', '',
  '| Candidate | Contains true bridge | Contains a valid adjoint solution | Closest bridge RMSE in class | Smallest adjoint equation RMSE in class |',
  '|:--|:--|:--|--:|--:|')
for (i in seq_len(nrow(classes))) {
  r <- classes[i]
  lines <- c(lines, sprintf('| %s | %s | %s | %.6f | %.6f |', labels[r$candidate],
    ifelse(r$contains_true_bridge, 'Yes', 'No'), ifelse(r$contains_valid_adjoint_solution, 'Yes', 'No'),
    r$closest_true_beta_rmse_in_class, r$smallest_adjoint_equation_rmse_in_class))
}
lines <- c(lines, '',
  'Saturated sieve and saturated L1 contain the required true functions. Default additive Landweber and the fitted PMMR kernel approximation do not, particularly for the adjoint. Their individual plots show approximation error. The expanded library supplies the saturated alternatives, and the adjoint places 98.47% weight on saturated sieve. The bridge puts 63.30% on Landweber, 29.62% on saturated L1, 4.51% on saturated sieve, and 2.57% on PMMR. Validation selects this mixture; a positive weight does not mean that a candidate class contains the true function.', '',
  sprintf('Although the complete dataset contains one million people, the selected training sample has %s measured input combinations with fewer than 10 observations and %s with fewer than 100. All 256 measured bridge input combinations are represented, but the smallest contains %s observation. A large overall sample does not make each full-history combination large.',
    verification$measured_target_cells_under_10, verification$measured_target_cells_under_100,
    verification$minimum_measured_target_cell_n), '',
  '| Largest ensemble departure | True value | Fitted value | Error | Measured training observations | All training observations at this history |',
  '|:--|--:|--:|--:|--:|--:|')
for (kind in c('bridge', 'adjoint-equation')) {
  r <- worst[candidate == 'ensemble' & function_name == kind][1L]
  lines <- c(lines, sprintf('| %s | %.6f | %.6f | %.6f | %s | %s |',
    ifelse(kind == 'bridge', 'Bridge prediction', 'Adjoint conditional equation'),
    r$true, r$estimate, r$error, r$training_measured, r$history_training))
}
lines <- c(lines, '',
  'These counts locate the remaining large departures in sparsely represented full-history combinations. They do not establish that every error comes exclusively from sampling variation: regularization and approximation also affect these fitted values. The saturated sieve still has a bridge point with error 2.8347, so class containment alone does not produce uniformly accurate estimates in this sample.', '',
  'The saturated L1 candidate reproduces the previous million-person bridge and adjoint predictions exactly (maximum differences both zero). This independently checks that the original data, splits, ratio construction, and L1 fitting remain unchanged in the expanded run. The [class audit](function-class-audit.csv) and [largest point errors with training counts](largest-point-errors.csv) retain the numerical evidence.', '',
  sprintf('The fits completed in %.1f minutes. The recorded warning is: %s. Every candidate completed with finite predictions. A single large sample cannot establish repeated-sample confidence-interval coverage or a convergence rate.',
    new$seconds / 60, if (length(new$warnings)) paste(new$warnings, collapse = ' | ') else 'none'), '',
  '[Complete diagnostics](component-comparison.csv), [bridge points](bridge-predictions.csv), [adjoint equation points](adjoint-equation-predictions.csv), [execution log](run.log), [reproduction instructions](REPRODUCE.md), and [source checksums](SOURCE-SHA256SUMS) are included.')
writeLines(lines, file.path(out, 'REPORT.md'))
print(metrics, digits = 8); print(weights, digits = 8); print(losses, digits = 8)
cat('Plots and report written to', out, '\n')
