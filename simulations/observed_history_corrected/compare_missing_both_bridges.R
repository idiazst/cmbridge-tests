# Read-only comparison of saved/recovered bridge fits and their scoring rule.
# Exact functions are diagnostic references, never fitting inputs.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
library(ggplot2)
root <- Sys.getenv('SIM_OUTPUT')
out <- file.path(root, 'diagnostics/binary-numerical-bridge-comparison')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
comparisons <- scores <- counts_out <- points_out <- list()
for (mechanism in c('binary_longitudinal', 'discrete_dose')) {
  tag <- if (mechanism == 'binary_longitudinal') 'binary-dose' else 'discrete-dose'
  label <- if (mechanism == 'binary_longitudinal') 'Binary treatment' else 'Numerical dose'
  x <- readRDS(file.path(root, sprintf('jobs/%s-n4000-r006.rds', mechanism)))
  saved <- readRDS(file.path(root, sprintf('diagnostics/%s-n4000-r006-sieve/recovered-sieve-fits.rds', tag)))
  g <- make_mechanism(mechanism, dose_max = max(2L, x$design$dose_max), visits = x$design$visits)
  p <- x$population
  stopifnot(identical(p, enumerate_data(g)))
  pt <- study_task(p, g)$task
  task <- study_task(saved$data, g, folds = x$folds, learner_folds = x$learner_folds,
    learner_groups = x$design$learner_groups, balance_measurement = x$design$balance_measurement)$task
  H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]
  B <- cell_key(pt$natural[, c(H, A), drop = FALSE])
  sample_B <- cell_key(task$natural[, c(H, A), drop = FALSE])
  Hkey <- cell_key(pt$natural[, H, drop = FALSE])
  measured <- p$R3 == 1L
  truth <- population_truth_functions(p, g)$beta[, 2L]
  # A diagnostic reference constant in C3 within each H2. It uses exact
  # observation probabilities only to examine the loss, not to fit a model.
  flat <- 1 / unname(conditional_mean(as.numeric(measured), p$probability, Hkey)[Hkey])
  references <- list('True bridge' = truth, 'Bridge constant in C3 within H2' = flat)
  moments <- lapply(names(references), function(reference) {
    residual <- rep(-1, nrow(p))
    residual[measured] <- references[[reference]][measured] - 1
    data.table(B = B, probability = p$probability, residual = residual)[, {
      mean <- sum(probability * residual) / sum(probability)
      .(mass = sum(probability), mean = mean,
        variance = max(0, sum(probability * residual^2) / sum(probability) - mean^2))
    }, by = B][, reference := reference]
  })
  moments <- rbindlist(moments)
  stopifnot(max(abs(moments[reference == 'True bridge']$mean)) < 1e-12)
  for (j in seq_along(x$folds)) {
    for (method in c('bridge', 'sieve')) {
      z <- fread(file.path(root, sprintf('diagnostics/%s-n4000-r006-%s/bridge-predictions.csv', tag, method)))[fold == j]
      z[, `:=`(mechanism = label, method = if (method == 'bridge') 'Ensemble' else 'Saturated sieve_md')]
      points_out[[length(points_out) + 1L]] <- z
      for (r in c(0L, 1L)) {
        zz <- z[R2 == r]
        comparisons[[length(comparisons) + 1L]] <- data.table(mechanism = label,
          method = z$method[1L], fold = j, R2 = r,
          rmse = sqrt(weighted.mean((zz$estimate - zz$true)^2, zz$probability)),
          signed_error = weighted.mean(zz$estimate - zz$true, zz$probability))
      }
    }
    rows <- x$folds[[j]]$training_set
    labels <- x$learner_folds[[j]]
    for (validation_label in sort(unique(labels))) {
      ii <- rows[labels == validation_label]
      cell_counts <- data.table(B = sample_B[ii])[, .(count = .N), by = B]
      counts_out[[length(counts_out) + 1L]] <- data.table(mechanism = label, fold = j,
        learner_fold = validation_label, validation_n = length(ii),
        represented_equation_cells = nrow(cell_counts),
        singleton_cells = sum(cell_counts$count == 1L),
        median_cell_count = median(cell_counts$count))
      for (reference_value in names(references)) {
        values <- merge(cell_counts, moments[reference == reference_value], by = 'B')
        stopifnot(nrow(values) == nrow(cell_counts))
        mean_term <- sum(values$count * values$mean^2) / length(ii)
        variance_term <- sum(values$variance) / length(ii)
        scores[[length(scores) + 1L]] <- data.table(mechanism = label, fold = j,
          learner_fold = validation_label, reference = reference_value,
          conditional_mean_squared_term = mean_term,
          residual_variance_term = variance_term,
          expected_score_independent_rows_given_cell_counts = mean_term + variance_term,
          max_population_bridge_equation_error = max(abs(moments[reference == reference_value]$mean)))
      }
    }
  }
}
points <- rbindlist(points_out)
binary <- points[mechanism == 'Binary treatment' & R2 == 1L]
binary[, `:=`(method = factor(method, levels = c('Ensemble', 'Saturated sieve_md')),
  fit = factor(paste('Training fit', fold), levels = paste('Training fit', 1:3)))]
figure <- ggplot(binary, aes(true, estimate)) +
  geom_abline(slope = 1, intercept = 0, color = '#3265a8', linewidth = .75) +
  geom_point(aes(size = probability), color = '#254a40', alpha = .65) +
  scale_size_area(max_size = 5, guide = 'none') +
  coord_equal(xlim = c(.8, 6.2), ylim = c(.8, 6.2), expand = FALSE) +
  facet_grid(method ~ fit) +
  labs(title = 'Binary treatment: true versus estimated bridge for outcome at time 3',
    subtitle = 'n = 4,000; replication 6; seed 5103006. R₂ = 1: the bridge solution is unique.',
    x = expression('True '*beta[2](H[2], C[3])),
    y = expression('Estimated '*beta[2](H[2], C[3])),
    caption = 'Blue line: y = x. Point area: population probability among measured final outcomes.') +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = 'bold'),
    plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0),
    plot.margin = margin(12, 16, 12, 12))
ggsave(file.path(out, 'binary-ensemble-and-sieve.png'), figure, width = 11.5, height = 8.3, dpi = 190, bg = 'white')
pdf_device <- if (capabilities('aqua')) {
  function(filename, ...) grDevices::quartz(type = 'pdf', file = filename, ...)
} else grDevices::cairo_pdf
ggsave(file.path(out, 'binary-ensemble-and-sieve.pdf'), figure,
  width = 11.5, height = 8.3, device = pdf_device, bg = 'white')
fwrite(rbindlist(comparisons), file.path(out, 'bridge-error-comparison.csv'))
fwrite(rbindlist(scores), file.path(out, 'validation-criterion-expectation.csv'))
fwrite(rbindlist(counts_out), file.path(out, 'validation-cell-counts.csv'))
fwrite(points, file.path(out, 'plotted-values.csv'))
cat('R2=1 bridge errors\n')
print(rbindlist(comparisons)[R2 == 1L])
cat('Expected loss for independent validation rows with the recorded cell counts\n')
print(rbindlist(scores)[, .(conditional_mean_squared_term = mean(conditional_mean_squared_term),
  residual_variance_term = mean(residual_variance_term),
  expected_score = mean(expected_score_independent_rows_given_cell_counts)), by = .(mechanism, reference)])
cat('These expectations audit the loss, not the original candidate CV scores. No estimator was fitted.\n')
