# Recreate only the deterministic sieve bridge candidate on the original data.
study_source <- Sys.getenv('STUDY_SOURCE', 'simulations/observed_history_corrected')
source(file.path(study_source, 'study.R'))
library(ggplot2)
root <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/missing-both-study-v2')
mechanism <- Sys.getenv('DIAGNOSTIC_MECHANISM', 'discrete_dose')
replicate <- as.integer(Sys.getenv('DIAGNOSTIC_REPLICATE', '6'))
label <- if (mechanism == 'binary_longitudinal') 'Binary treatment' else 'Numerical dose'
tag <- if (mechanism == 'binary_longitudinal') 'binary-dose' else 'discrete-dose'
out <- file.path(root, sprintf('diagnostics/%s-n4000-r%03d-sieve', tag, replicate))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(file.path(root, sprintf('jobs/%s-n4000-r%03d.rds', mechanism, replicate)))
stopifnot(x$design$n == 4000L, x$design$replicate == replicate,
  identical(x$design$mechanism, mechanism), identical(x$design$visits, 'missing_both'))
set.seed(x$design$seed)
g <- make_mechanism(mechanism, dose_max = max(2L, x$design$dose_max), visits = x$design$visits)
d <- draw_data(x$design$n, g)
prepared <- study_task(d, g, folds = x$folds, learner_folds = x$learner_folds,
  learner_groups = x$design$learner_groups, balance_measurement = x$design$balance_measurement)
stopifnot(identical(prepared$task$folds, x$folds),
  identical(prepared$task$learner_folds, x$learner_folds))
task <- prepared$task
args <- study_arguments(g)
# Check regeneration against the saved counts for both outcomes and all folds.
for (j in seq_along(x$folds)) for (s in 1:2) {
  train <- x$folds[[j]]$training_set
  H <- task$vars$history('A', s); A <- task$vars$A[[s]]
  B <- cell_key(task$natural[train, c(H, A), drop = FALSE])
  V <- cell_key(cbind(task$natural[train, H, drop = FALSE],
    d[train, args$health[[s]], drop = FALSE]))
  M <- d[[args$measurement[s]]][train]
  for (kind_value in c('B', 'V')) {
    keys <- if (kind_value == 'B') B else V
    got <- data.table(cell = keys, M = M)[,
      .(training = .N, measured = sum(M)), by = cell][order(cell)]
    original <- x$occupancy[fold == j & horizon == s & kind == kind_value,
      .(cell, training, measured)][order(cell)]
    stopifnot(isTRUE(all.equal(got, original, check.attributes = FALSE)))
  }
}
p <- x$population
stopifnot(identical(p, enumerate_data(g)))
truth <- population_truth_functions(p, g)
pt <- study_task(p, g)$task
H <- task$vars$history('A', 2L); A <- task$vars$A[[2L]]
B <- as.matrix(task$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(task$natural[, H, drop = FALSE], d[, args$health[[2L]], drop = FALSE]))
observed <- p$R3 == 1L
newV <- as.matrix(cbind(pt$natural[observed, H, drop = FALSE],
  p[observed, args$health[[2L]], drop = FALSE]))
key <- cell_key(newV)
candidate <- cm_library()[['sieve_md']]
stopifnot(candidate$method == 'sieve_md', candidate$control$target_basis == 'cell',
  candidate$control$instrument_basis == 'cell')
models <- points <- metrics <- list()
for (j in seq_along(x$folds)) {
  train <- x$folds[[j]]$training_set
  model <- cmbridge::fit_bridge(B[train, , drop = FALSE], V[train, , drop = FALSE],
    d$R3[train], method = candidate$method, control = candidate$control)
  estimate <- predict(model, newV)
  z <- data.table(cell = key, R2 = p$R2[observed],
    true = truth$beta[observed, 2L], estimate = estimate,
    probability = p$probability[observed])[, {
      stopifnot(diff(range(true)) < 1e-12, diff(range(estimate)) < 1e-12)
      .(R2 = R2[1L], true = true[1L], estimate = estimate[1L],
        probability = sum(probability))
    }, by = cell]
  z[, fold := j]
  stopifnot(all(is.finite(z$estimate)))
  residual <- rep(-1, nrow(p)); residual[observed] <- estimate - 1
  moments <- conditional_mean(residual, p$probability,
    cell_key(pt$natural[, c(H, A), drop = FALSE]))
  ensemble <- x$fitted_functions[[paste('sdr', 'beta1_lambda1_ratio1_m1', j, sep = '/')]]$beta[observed, 2L]
  metrics[[j]] <- data.table(fold = j, training_n = length(train),
    rmse_sieve = sqrt(weighted.mean((estimate - truth$beta[observed, 2L])^2, p$probability[observed])),
    rmse_ensemble = sqrt(weighted.mean((ensemble - truth$beta[observed, 2L])^2, p$probability[observed])),
    bridge_equation_max_error = max(abs(moments)),
    lower_bound_fraction = weighted.mean(estimate <= 1 + 1e-6, p$probability[observed]),
    upper_bound_fraction = weighted.mean(estimate >= 6 - 1e-6, p$probability[observed]))
  models[[j]] <- model; points[[j]] <- z
}
points <- rbindlist(points)
points[, `:=`(fit_label = factor(paste('Training fit', fold), levels = paste('Training fit', 1:3)),
  measurement_label = factor(ifelse(R2 == 1L, 'R₂ = 1', 'R₂ = 0'), levels = c('R₂ = 1', 'R₂ = 0')))]
limits <- range(points$true, points$estimate)
limits <- limits + c(-1, 1) * max(.025, diff(limits) * .04)
figure <- ggplot(points, aes(true, estimate)) +
  geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
  geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
  scale_size_area(max_size = 5, guide = 'none') +
  coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
  facet_grid(measurement_label ~ fit_label) +
  labs(title = 'Saturated sieve_md alone: true versus estimated bridge',
    subtitle = sprintf('%s; n = 4,000; replication %d; seed %d. Original training splits.',
      label, replicate, x$design$seed),
    x = expression('True '*beta[2](H[2], C[3])),
    y = expression('Estimated '*beta[2](H[2], C[3])),
    caption = paste('Blue line: y = x. Point area: population probability among measured outcomes.',
      'R₂ = 0: bridge solutions are not unique; the x-axis uses the selected valid DGP solution.', sep = '\n')) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
    plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 10),
    plot.margin = margin(12, 16, 12, 12))
ggsave(file.path(out, 'sieve-true-versus-estimated.png'), figure,
  width = 11.5, height = 8.5, dpi = 190, bg = 'white')
pdf_device <- if (capabilities('aqua')) {
  function(filename, ...) grDevices::quartz(type = 'pdf', file = filename, ...)
} else grDevices::cairo_pdf
ggsave(file.path(out, 'sieve-true-versus-estimated.pdf'), figure,
  width = 11.5, height = 8.5, device = pdf_device, bg = 'white')
fwrite(points, file.path(out, 'bridge-predictions.csv'))
fwrite(rbindlist(metrics), file.path(out, 'fit-diagnostics.csv'))
saveRDS(list(models = models, folds = x$folds, learner_folds = x$learner_folds,
  data = d, design = x$design, candidate = candidate), file.path(out, 'recovered-sieve-fits.rds'))
writeLines(c('# Sieve-only bridge for the selected dataset', '',
  'The original ensemble checkpoint did not retain individual candidate predictions. This check reconstructs only its deterministic saturated sieve_md bridge candidate using the frozen package and exactly the original seeded dataset, predictors, controls and outer training splits. The regenerated data agree with the saved predictor-combination counts for both outcomes and every fold.', '',
  'No treatment ratios, adjoints, sequential regressions or final SDR/TMLE estimators were refitted. The full simulation remains stopped.', '',
  'The plotted function is the bridge for outcome time 3. The figure uses the same DGP reference as the ensemble plot. At R₂ = 0 the solution is not unique, so differences from this reference alone do not establish an equation error.', '',
  '![Sieve-only bridge](sieve-true-versus-estimated.png)'), file.path(out, 'README.md'))
cat('Recreated three sieve bridge fits on the original dataset and splits.\n')
print(rbindlist(metrics))
