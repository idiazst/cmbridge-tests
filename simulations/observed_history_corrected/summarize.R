library(data.table)
library(ggplot2)
outdir <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/missing-both-study-v2")
paths <- list.files(file.path(outdir, "jobs"), "\\.rds$", full.names = TRUE)
if (!length(paths)) stop("No completed checkpoints.")
design <- fread(file.path(outdir, "design.csv"))
results <- diagnostics <- errors <- tuning <- weights <- logs <- list()
candidate_failures <- penalty_trial_failures <- list()
for (path in paths) {
  x <- readRDS(path); d <- as.data.table(x$design)
  if (!is.null(x$results) && nrow(x$results)) {
    if (length(x$job_warnings)) x$results[, warnings := paste(warnings, paste(x$job_warnings, collapse = " | "), sep = " | ")]
    results[[length(results) + 1L]] <- x$results
  }
  if (!is.null(x$diagnostics) && nrow(x$diagnostics)) diagnostics[[length(diagnostics) + 1L]] <- x$diagnostics
  if (!is.null(x$errors) && nrow(x$errors)) errors[[length(errors) + 1L]] <- cbind(d, x$errors)
  if (!is.null(x$tuning) && nrow(x$tuning)) tuning[[length(tuning) + 1L]] <- cbind(d, x$tuning)
  if (!is.null(x$ensemble_weights) && nrow(x$ensemble_weights))
    weights[[length(weights) + 1L]] <- cbind(d, x$ensemble_weights)
  if (!is.null(x$candidate_failures) && nrow(x$candidate_failures))
    candidate_failures[[length(candidate_failures) + 1L]] <- x$candidate_failures
  if (!is.null(x$penalty_trial_failures) && nrow(x$penalty_trial_failures))
    penalty_trial_failures[[length(penalty_trial_failures) + 1L]] <- x$penalty_trial_failures
  logs[[length(logs) + 1L]] <- cbind(d, data.table(seconds = x$seconds,
    comparisons_successful = if (is.null(x$results) || !nrow(x$results)) 0L else uniqueN(x$results[, .(scenario, estimator)]),
    errors = if (is.null(x$errors)) 0L else nrow(x$errors),
    warnings = paste(x$job_warnings, collapse = " | ")))
}
log <- rbindlist(logs, fill = TRUE)
fwrite(log, file.path(outdir, "job_log.csv"))
for (item in list(list(candidate_failures, 'nested-candidate-failures.csv'),
                  list(penalty_trial_failures, 'nested-penalty-trial-failures.csv'))) {
  recorded <- rbindlist(item[[1L]], fill = TRUE)
  if (!ncol(recorded)) recorded <- data.table(error = character())
  fwrite(recorded, file.path(outdir, item[[2L]]))
}
timing <- log[, .(datasets_timed = .N, mean_seconds = mean(seconds), median_seconds = median(seconds)),
              by = .(mechanism, n)]
timing <- merge(design[, .(planned = .N), by = .(mechanism, n)], timing,
                by = c("mechanism", "n"), all.x = TRUE)
timing[, remaining_datasets := planned - ifelse(is.na(datasets_timed), 0L, datasets_timed)]
timing[, estimated_remaining_worker_hours := remaining_datasets * mean_seconds / 3600]
fwrite(timing, file.path(outdir, "execution-time-estimate.csv"))
error_table <- rbindlist(errors, fill = TRUE)
if (!ncol(error_table)) error_table <- data.table(error = character())
fwrite(error_table, file.path(outdir, "errors.csv"))
x <- rbindlist(results, fill = TRUE)
if (!nrow(x)) stop("Checkpoints contain no successful comparisons; inspect errors.csv.")
x[, estimation_error := estimate - truth]
x[, covered := lower <= truth & upper >= truth]
fwrite(x, file.path(outdir, "replicates.csv"))
keys <- c("mechanism", "horizon", "estimator", "scenario", "beta_correct", "lambda_correct",
          "ratio_correct", "m_correct", "sufficient_consistency", "n")
component_columns <- c('mean_sequential_eif','mean_bridge_adjoint_eif','mean_total_eif')
if(all(component_columns %in% names(x))) {
  stopifnot(max(abs(x$mean_sequential_eif+x$mean_bridge_adjoint_eif-x$mean_total_eif))<1e-10,
    max(abs(x$mean_total_eif-x$estimation_error))<1e-10)
  fwrite(x[,c(keys,'replicate','seed',component_columns),with=FALSE],
    file.path(outdir,'eif-component-means-by-replication.csv'))
  component_summary <- x[,.(successful=.N,
    mean_sequential_eif=mean(mean_sequential_eif),
    mean_bridge_adjoint_eif=mean(mean_bridge_adjoint_eif),
    mean_total_eif=mean(mean_total_eif)),by=keys]
  fwrite(component_summary,file.path(outdir,'eif-component-means.csv'))
}
summary <- x[, .(successful = .N, truth = mean(truth), mean_estimate = mean(estimate),
  bias = mean(estimation_error), bias_mcse = sd(estimation_error) / sqrt(.N),
  root_n_bias = sqrt(n[1]) * mean(estimation_error),
  root_n_bias_mcse = sqrt(n[1]) * sd(estimation_error) / sqrt(.N),
  empirical_sd = sd(estimate), rmse = sqrt(mean(estimation_error^2)),
  mean_se = mean(se), mean_interval_length = mean(upper - lower),
  coverage = mean(covered), coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
  joint_curve_coverage = mean(joint_curve_covered),
  replications_with_warnings = sum(nzchar(warnings))), by = keys]
planned <- design[, .(planned = .N), by = .(mechanism, n)]
completed <- log[, .(datasets_completed = .N), by = .(mechanism, n)]
spec_file <- file.path(outdir, "specifications.csv")
specifications <- if (file.exists(spec_file)) fread(spec_file) else
  CJ(estimator = c("sdr", "tmle"), beta_correct = c(TRUE, FALSE),
    lambda_correct = c(TRUE, FALSE), ratio_correct = c(TRUE, FALSE), m_correct = c(TRUE, FALSE))
specifications[, join := 1L]; planned[, join := 1L]
expected <- merge(planned, specifications, by = "join", allow.cartesian = TRUE)
expected[, join := NULL]; planned[, join := NULL]
expected <- expected[, .(horizon = seq_len(if (mechanism == "binary_point") 1L else 2L)),
  by = setdiff(names(expected), "horizon")]
expected[, scenario := paste0("beta", as.integer(beta_correct), "_lambda", as.integer(lambda_correct),
  "_ratio", as.integer(ratio_correct), "_m", as.integer(m_correct))]
expected[, sufficient_consistency := (beta_correct | lambda_correct) & (ratio_correct | m_correct)]
summary <- merge(expected, summary, by = keys, all.x = TRUE)
summary <- merge(summary, completed, by = c("mechanism", "n"), all.x = TRUE)
summary[is.na(successful), successful := 0L]
summary[is.na(datasets_completed), datasets_completed := 0L]
summary[, `:=`(failed = datasets_completed - successful, pending = planned - datasets_completed,
               se_calibration = mean_se / empirical_sd)]
z <- qnorm(.975)
summary[, coverage_mc_lower := (coverage + z^2 / (2 * successful) -
  z * sqrt(coverage * (1 - coverage) / successful + z^2 / (4 * successful^2))) / (1 + z^2 / successful)]
summary[, coverage_mc_upper := (coverage + z^2 / (2 * successful) +
  z * sqrt(coverage * (1 - coverage) / successful + z^2 / (4 * successful^2))) / (1 + z^2 / successful)]
fwrite(summary, file.path(outdir, "summary.csv"))
dg <- rbindlist(diagnostics, fill = TRUE)
fwrite(dg, file.path(outdir, "population_diagnostics.csv"))
if (nrow(dg)) {
  metrics <- c("population_bias", "sequential_remainder", "bridge_remainder", "sequential_product_norm",
    "beta_mse_measured", "lambda_mse_to_selected_solution", "ratio_mse_sum", "regression_mse_sum",
    "bridge_equation_max_error", "adjoint_equation_max_error", "beta_lower_fraction", "beta_upper_fraction")
  if (!"fold_weight" %in% names(dg)) dg[, fold_weight := 1 / .N,
    by = .(mechanism, n, replicate, estimator, scenario, horizon)]
  dd <- dg[, lapply(.SD, function(value) sum(value * fold_weight)),
    by = .(mechanism, n, replicate, estimator, scenario, horizon), .SDcols = metrics]
  ds <- melt(dd, id.vars = c("mechanism", "n", "replicate", "estimator", "scenario", "horizon"))[
    , .(mean = mean(value), sd = sd(value), q025 = quantile(value, .025), q975 = quantile(value, .975)),
    by = .(mechanism, n, estimator, scenario, horizon, variable)]
  fwrite(ds, file.path(outdir, "population_diagnostic_summary.csv"))
}
tt <- rbindlist(tuning, fill = TRUE)
if (nrow(tt)) {
  fwrite(tt, file.path(outdir, "penalty_and_weight_selections.csv"))
  fwrite(tt[selected == TRUE, .(selected_fits = .N),
    by = .(mechanism, n, kind, scale)], file.path(outdir, "penalty_selection_summary.csv"))
}
ww <- rbindlist(weights, fill = TRUE)
if (nrow(ww)) {
  fwrite(ww, file.path(outdir, "ensemble_weights.csv"))
  weight_keys <- intersect(c("mechanism", "n", "kind", "horizon_or_depth", "candidate", "method", "estimator"), names(ww))
  ws <- ww[, .(fits = .N, mean_weight = mean(weight), median_weight = median(weight),
               positive_weight_fraction = mean(weight > 1e-8)), by = weight_keys]
  fwrite(ws, file.path(outdir, "ensemble_weight_summary.csv"))
}
figures <- file.path(outdir, "figures")
dir.create(figures, showWarnings = FALSE)
save_plot <- function(p, name, width = 12, height = 8) {
  p <- p + theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(figures, paste0(name, ".png")), p, width = width, height = height, dpi = 170)
  pdf_device <- if (capabilities("aqua")) {
    function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...)
  } else grDevices::cairo_pdf
  ggsave(file.path(figures, paste0(name, ".pdf")), p, width = width, height = height,
         device = pdf_device)
}
# The final outcome is the primary comparison within each mechanism.
last <- summary[horizon == ifelse(mechanism == "binary_point", 1, 2)]
last[, specification := paste0("β ", ifelse(beta_correct, "correct", "constant"), "; λ ",
  ifelse(lambda_correct, "correct", "constant"), "\nω ", ifelse(ratio_correct, "correct", "omits H"),
  "; m ", ifelse(m_correct, "correct", "constant"))]
for (name in unique(last$mechanism)) {
  y <- last[mechanism == name]
  panels <- uniqueN(y$specification); columns <- if (panels == 1L) 1L else 4L
  width <- if (panels == 1L) 7 else 14; height <- if (panels == 1L) 5 else 11
  base <- ggplot(y, aes(n, color = estimator, group = estimator)) + scale_x_log10(breaks = sort(unique(design$n)))
  save_plot(base + aes(y = bias) + geom_hline(yintercept = 0, color = "grey60") + geom_line() + geom_point() +
    geom_errorbar(aes(ymin = bias - 1.96 * bias_mcse, ymax = bias + 1.96 * bias_mcse), width = 0, na.rm = TRUE) +
    facet_wrap(~specification, ncol = columns, scales = "free_y") + labs(title = paste(name, "— estimation error"),
    subtitle = "Bars show simulation uncertainty in the mean error.", x = "Sample size", y = "Mean estimate − truth"), paste0(name, "-bias"), width, height)
  save_plot(base + aes(y = root_n_bias) + geom_hline(yintercept = 0, color = "grey60") + geom_line() + geom_point() +
    facet_wrap(~specification, ncol = columns, scales = "free_y") + labs(title = paste(name, "— √n × mean estimation error"),
    x = "Sample size", y = "√n × bias"), paste0(name, "-root-n-bias"), width, height)
  save_plot(base + aes(y = coverage) + geom_hline(yintercept = .95, linetype = 2) + geom_line() + geom_point() +
    geom_errorbar(aes(ymin = coverage_mc_lower, ymax = coverage_mc_upper), width = 0) +
    facet_wrap(~specification, ncol = columns) + coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste(name, "— coverage of 95% intervals"), x = "Sample size", y = "Coverage"), paste0(name, "-coverage"), width, height)
  save_plot(ggplot(y, aes(empirical_sd, mean_se, color = estimator, shape = factor(n))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(na.rm = TRUE) + facet_wrap(~specification, ncol = columns, scales = "free") +
    labs(title = paste(name, "— standard-error calibration"), x = "SD across replications", y = "Mean estimated SE", shape = "n"),
    paste0(name, "-standard-errors"), width, height)
}
correct <- summary[beta_correct & lambda_correct & ratio_correct & m_correct]
save_plot(ggplot(correct, aes(horizon + 1, mean_estimate, color = estimator, group = estimator)) +
  geom_line() + geom_point() + geom_line(aes(y = truth), color = "black", linetype = 2) +
  facet_grid(mechanism ~ n) + scale_x_continuous(breaks = 2:3) +
  labs(title = "Estimated curves with all four function classes correctly specified", x = "Outcome time", y = "Mean estimate; dashed line is truth"),
  "correctly-specified-curves", 15, 7)
complete <- nrow(log) == nrow(design)
status <- c(if(file.exists(file.path(outdir,'run-stopped-at.txt'))) 'Stopped at the user\'s request.' else NULL,
  paste("Completed datasets:", nrow(log), "of", nrow(design)),
  paste("Successful comparisons:", nrow(unique(x[, .(mechanism, n, replicate, estimator, scenario)]))),
  paste("Datasets with reported errors:", sum(log$errors > 0)),
  paste("Remaining datasets:", nrow(design) - nrow(log)),
  if (complete) "All planned datasets processed. Successful fits, failures, and warnings are retained." else
  "Coverage and bias conclusions require the planned replication counts; partial outputs are preliminary.")
writeLines(status, file.path(outdir, "STATUS.txt")); cat(status, sep = "\n")
