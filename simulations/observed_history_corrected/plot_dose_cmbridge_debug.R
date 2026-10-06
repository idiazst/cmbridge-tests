# Redraw completed component checks with common x/y limits for y = x.
# This reads saved predictions only; it does not fit any learner.
library(data.table)
library(ggplot2)
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
bp <- fread(file.path(out, 'bridge-predictions.csv'))
ap <- fread(file.path(out, 'adjoint-equation-predictions.csv'))
summary <- fread(file.path(out, 'fit-diagnostics.csv'))
stopifnot(all(is.finite(bp$true)), all(is.finite(bp$estimate)),
          all(is.finite(ap$true)), all(is.finite(ap$estimate)),
          all(bp$probability > 0), all(ap$probability > 0))
reference <- bp[n == 4000L & seed == 5203002L & fold == 2L]
large <- bp[source == 'large check']
large_a <- ap[source == 'large check']
stopifnot(nrow(reference) == 256L, nrow(large) == 256L, nrow(large_a) == 256L)

limits <- function(points) {
  r <- range(points$true, points$estimate)
  r + c(-1, 1) * max(diff(r) * .04, .025)
}
bridge_limits <- limits(rbind(reference, large))
base_plot <- function(points, title, subtitle, axis_limits,
                      xlab = 'True value from the DGP', ylab = 'Estimated value') {
  ggplot(points, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
    geom_point(aes(size = probability), alpha = .60, color = '#254a40') +
    scale_size_area(max_size = 5, guide = 'none') +
    coord_equal(xlim = axis_limits, ylim = axis_limits, expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         caption = 'Blue line: y = x. Point area: generating probability among measured outcomes.') +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(size = 10.5, margin = margin(b = 12)),
          plot.caption = element_text(size = 9, hjust = 0),
          strip.text = element_text(size = 11, face = 'bold'),
          plot.margin = margin(12, 14, 12, 12))
}
save_plot <- function(plot, name, width = 7, height = 6.4) {
  for (extension in c('png', 'pdf')) {
    ggsave(file.path(out, 'figures', paste0(name, '.', extension)), plot,
           width = width, height = height, dpi = 200, bg = 'white')
  }
}
save_plot(base_plot(reference, 'Bridge at time 2: estimate versus truth',
  'n = 4,000; seed 5203002; training sample 2 (2,667 people)', bridge_limits),
  'bridge-y-equals-x-reference')
large_summary <- summary[source == 'large check']
counts <- fread(file.path(out, 'large-training-counts.csv'))
save_plot(base_plot(large, 'Bridge at time 2: estimate versus truth',
  sprintf('n = %s; seed %s; training sample 2 (%s people)',
          format(large_summary$n, big.mark = ','), large_summary$seed,
          format(sum(counts[time == 2L]$training), big.mark = ',')), bridge_limits),
  'bridge-y-equals-x-large')

comparison <- rbind(reference, large)
comparison[, label := ifelse(source == 'large check', 'n = 1,000,000', 'n = 4,000')]
comparison[, label := factor(label, levels = c('n = 4,000', 'n = 1,000,000'))]
save_plot(base_plot(comparison, 'Bridge at time 2: estimate versus truth',
  'Same DGP and learner; positive L1 penalty selected by cross-validation.', bridge_limits) +
  facet_wrap(~label, ncol = 2), 'bridge-y-equals-x-comparison', 11, 6.4)

save_plot(base_plot(large_a, 'Adjoint at time 2: conditional equation',
  'n = 1,000,000; assess the equation because the adjoint is not unique.', limits(large_a),
  'True conditional mean from the DGP', 'Conditional mean of fitted adjoint'),
  'adjoint-equation-y-equals-x-large')
reference_a <- ap[n == 4000L & seed == 5203002L & fold == 2L]
comparison_a <- rbind(reference_a, large_a)
comparison_a[, label := ifelse(source == 'large check', 'n = 1,000,000', 'n = 4,000')]
comparison_a[, label := factor(label, levels = c('n = 4,000', 'n = 1,000,000'))]
save_plot(base_plot(comparison_a, 'Adjoint at time 2: conditional equation',
  'Same DGP and learner; treatment ratios estimated within the training samples.', limits(comparison_a),
  'True conditional mean from the DGP', 'Conditional mean of fitted adjoint') +
  facet_wrap(~label, ncol = 2), 'adjoint-equation-y-equals-x-comparison', 11, 6.4)

# Confirm that the plotted values reproduce the table's probability-weighted errors.
for (points in list(reference, large)) {
  row <- summary[n == points$n[1] & seed == points$seed[1] &
                   fold == points$fold[1] & time == 2L]
  stopifnot(nrow(row) == 1L)
  rmse <- sqrt(sum(points$probability * (points$estimate - points$true)^2) /
                 sum(points$probability))
  stopifnot(abs(rmse - row$beta_rmse) < 1e-10)
}
for (points in list(reference_a, large_a)) {
  row <- summary[n == points$n[1] & seed == points$seed[1] &
                   fold == points$fold[1] & time == 2L]
  rmse <- sqrt(sum(points$probability * (points$estimate - points$true)^2) /
                 sum(points$probability))
  stopifnot(abs(rmse - row$adjoint_equation_rmse) < 1e-10)
}
cat('Figures redrawn; plotted values reproduce the numerical diagnostics.\n')
