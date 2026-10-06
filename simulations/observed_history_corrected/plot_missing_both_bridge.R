# Plot saved outer-training bridge predictions; no nuisance models are fitted.
study_source <- Sys.getenv('STUDY_SOURCE', 'simulations/observed_history_corrected')
source(file.path(study_source, 'study.R'))
library(ggplot2)
root <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/missing-both-study-v2')
mechanism <- Sys.getenv('DIAGNOSTIC_MECHANISM', 'discrete_dose')
replicate <- as.integer(Sys.getenv('DIAGNOSTIC_REPLICATE', '6'))
label <- if (mechanism == 'binary_longitudinal') 'Binary treatment' else 'Numerical dose'
tag <- if (mechanism == 'binary_longitudinal') 'binary-dose' else 'discrete-dose'
out <- file.path(root, sprintf('diagnostics/%s-n4000-r%03d-bridge', tag, replicate))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(root, sprintf('jobs/%s-n4000-r%03d.rds', mechanism, replicate))
x <- readRDS(checkpoint)
stopifnot(x$design$n == 4000L, x$design$replicate == replicate,
  identical(x$design$mechanism, mechanism), identical(x$design$visits, 'missing_both'))
g <- make_mechanism(mechanism, dose_max = max(2L, x$design$dose_max), visits = x$design$visits)
require_valid_mechanism(g)
p <- x$population
stopifnot(identical(p, enumerate_data(g)))
truth <- population_truth_functions(p, g)
pt <- study_task(p, g)$task
observed <- p$R3 == 1L
H <- pt$vars$history('A', 2L)
key <- cell_key(cbind(pt$natural[observed, H, drop = FALSE],
  p[observed, c('C3_covariate', 'Y3'), drop = FALSE]))
points <- metrics <- list()
for (j in seq_along(x$folds)) {
  path <- paste('sdr', 'beta1_lambda1_ratio1_m1', j, sep = '/')
  nu <- x$fitted_functions[[path]]
  stopifnot(identical(nu$beta,
    x$fitted_functions[[sub('^sdr/', 'tmle/', path)]]$beta))
  z <- data.table(cell = key, R2 = p$R2[observed],
    true = truth$beta[observed, 2L], estimate = nu$beta[observed, 2L],
    probability = p$probability[observed])
  z <- z[, {
    stopifnot(diff(range(true)) < 1e-12, diff(range(estimate)) < 1e-12,
      uniqueN(R2) == 1L)
    .(R2 = R2[1L], true = true[1L], estimate = estimate[1L],
      probability = sum(probability))
  }, by = cell]
  z[, fold := j]
  mse <- weighted.mean((z$estimate - z$true)^2, z$probability)
  reported <- x$diagnostics[estimator == 'sdr' & horizon == 2L & fold == j,
    beta_mse_measured]
  stopifnot(length(reported) == 1L, abs(mse - reported) < 1e-12,
    all(is.finite(z$estimate)), all(z$probability > 0))
  points[[j]] <- z
  metrics[[j]] <- data.table(fold = j,
    training_n = length(x$folds[[j]]$training_set),
    validation_n = length(x$folds[[j]]$validation_set),
    plotted_combinations = nrow(z), mse_measured = mse,
    max_absolute_error = max(abs(z$estimate - z$true)))
}
points <- rbindlist(points)
points[, `:=`(fit_label = factor(paste('Training fit', fold),
  levels = paste('Training fit', 1:3)),
  measurement_label = factor(ifelse(R2 == 1L, 'R₂ = 1', 'R₂ = 0'),
    levels = c('R₂ = 1', 'R₂ = 0')))]
limits <- range(points$true, points$estimate)
limits <- limits + c(-1, 1) * max(.025, diff(limits) * .04)
figure <- ggplot(points, aes(true, estimate)) +
  geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
  geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
  scale_size_area(max_size = 5, guide = 'none') +
  coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
  facet_grid(measurement_label ~ fit_label) +
  labs(title = 'Bridge for outcome at time 3: true versus estimated',
    subtitle = sprintf('%s; n = 4,000; replication %d; seed %d. Saved outer training fits.',
      label, replicate, x$design$seed),
    x = expression('True '*beta[2](H[2], C[3])),
    y = expression('Estimated '*beta[2](H[2], C[3])),
    caption = paste(
      'Blue line: y = x. Point area: population probability among measured outcomes.',
      'R₂ = 0: bridge solutions are not unique; the x-axis uses the selected valid DGP solution.',
      sep = '\n')) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
    plot.title = element_text(face = 'bold'),
    plot.caption = element_text(hjust = 0, size = 10),
    plot.margin = margin(12, 16, 12, 12))
ggsave(file.path(out, 'bridge-true-versus-estimated.png'), figure,
  width = 11.5, height = 8.5, dpi = 190, bg = 'white')
pdf_device <- if (capabilities('aqua')) {
  function(filename, ...) grDevices::quartz(type = 'pdf', file = filename, ...)
} else grDevices::cairo_pdf
ggsave(file.path(out, 'bridge-true-versus-estimated.pdf'), figure,
  width = 11.5, height = 8.5, device = pdf_device, bg = 'white')
fwrite(points, file.path(out, 'bridge-predictions.csv'))
fwrite(rbindlist(metrics), file.path(out, 'plot-validation.csv'))
fwrite(x$results[estimator == 'sdr' & horizon == 2L],
  file.path(out, 'selected-estimator-result.csv'))
writeLines(c(paste('# Bridge plot:', label), '',
  sprintf('This is replication %d, seed %d, n = 4,000, with missingness at both follow-ups.',
    replicate, x$design$seed), '',
  'The figure compares the saved bridge estimates for outcome time 3 with the selected valid DGP bridge, over the finite population support where the final outcome is measured. Each panel uses a separate original outer training fit. SDR and TMLE share these bridge fits. No fitting or new simulation was run.', '',
  'Combinations with the same bridge predictors are grouped; point areas show their generating probabilities among measured outcomes. The plotted mean squared errors reproduce the original saved diagnostics.', '',
  'At R₂ = 1 the conditional operator has full column rank and the bridge solution is unique. At R₂ = 0 multiple bridge solutions are valid, so distance from the selected reference alone is not evidence of an invalid solution.', '',
  '![True versus estimated bridge](bridge-true-versus-estimated.png)'),
  file.path(out, 'README.md'))
cat('Plotted saved bridge predictions; errors match checkpoint diagnostics.\n')
print(rbindlist(metrics))
