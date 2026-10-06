# Reconstruct predictions from saved positive-penalty coefficients and plot.
library(data.table)
library(ggplot2)
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
points <- fread(file.path(out, 'bridge-predictions.csv'))[source == 'large check']
record <- readRDS(file.path(out, 'positive-penalty-bridge-fits.rds'))
table <- fread(file.path(out, 'positive-penalty-shrinkage-comparison.csv'))
all <- list()
for (i in seq_len(nrow(record$scenarios))) {
  model <- record$models[[i]]
  index <- match(points$cell, model$target_levels)
  stopifnot(!anyNA(index))
  prediction <- pmax(model$tuning$lower, pmin(model$tuning$upper,
    model$coefficients[1L] + model$coefficients[1L + index]))
  got <- copy(points)
  set(got, j = 'estimate', value = prediction)
  got[, `:=`(scenario = record$scenarios$scenario[i], scale = record$scenarios$scale[i],
              tolerance = record$scenarios$tolerance[i])]
  rmse <- sqrt(weighted.mean((prediction - points$true)^2, points$probability))
  stopifnot(abs(rmse - table$beta_rmse[i]) < 1e-10)
  all[[i]] <- got
}
all <- rbindlist(all)
fwrite(all, file.path(out, 'positive-penalty-bridge-predictions.csv'))
plot_data <- all[tolerance == 1e-8]
plot_data[, label := factor(sprintf('Scale %s%s', format(scale),
  ifelse(scale == .01, ' (CV selected)', '')), levels = unique(sprintf('Scale %s%s', format(scale),
  ifelse(scale == .01, ' (CV selected)', ''))))]
limit <- range(plot_data$true, plot_data$estimate)
limit <- limit + c(-1, 1) * diff(limit) * .04
plot <- ggplot(plot_data, aes(true, estimate)) +
  geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .65) +
  geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
  scale_size_area(max_size = 4.2, guide = 'none') +
  coord_equal(xlim = limit, ylim = limit, expand = FALSE) + facet_wrap(~label, ncol = 3) +
  labs(title = 'Million-person bridge check: positive penalty comparison',
    subtitle = 'Same training sample and learner; scales are divided by the square root of 666,667 training people.',
    x = 'True value from the DGP', y = 'Estimated value',
    caption = 'Blue line: y = x. Point area: generating probability among measured outcomes. Truth is used only for evaluation.') +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'),
    plot.subtitle = element_text(size = 10, margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, hjust = 0),
    strip.text = element_text(face = 'bold'), plot.margin = margin(12,12,12,12))
for (extension in c('png','pdf')) ggsave(file.path(out,'figures',paste0('positive-penalty-comparison.',extension)),
  plot, width = 12, height = 5.5, dpi = 200, bg = 'white')

# Compare only the points far from the originally selected fit.
wide <- dcast(all[tolerance == 1e-8], cell + true + probability ~ scale, value.var = 'estimate')
original_error <- abs(wide[['0.01']] - wide$true)
selected <- wide[original_error > .5]
selected[, `:=`(selected_error = `0.01` - true, lower_penalty_error = `0.001` - true,
                lower_penalty_improves_absolute_error = abs(`0.001` - true) < abs(`0.01` - true))]
fwrite(selected[order(-abs(selected_error))], file.path(out, 'shrinkage-far-point-comparison.csv'))
cat('All plotted predictions reconstructed from fitted coefficients and matched the metric table.\n')
print(selected[, .(far_points = .N, improved_with_lower_penalty = sum(lower_penalty_improves_absolute_error))])
