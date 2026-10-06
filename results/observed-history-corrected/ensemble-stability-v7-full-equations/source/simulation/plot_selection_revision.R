library(data.table)
library(ggplot2)
root <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/selection-revision")
folders <- list.dirs(root, recursive = FALSE, full.names = TRUE)
folders <- folders[grepl("-n(4000|20000)-seed", basename(folders))]
read_tables <- function(pattern) rbindlist(lapply(folders, function(folder) {
  files <- list.files(folder, pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  d <- rbindlist(lapply(files, fread), fill = TRUE)
  d[, mechanism := sub("-n.*", "", basename(folder))]
  d[, n := as.integer(sub(".*-n([0-9]+)-.*", "\\1", basename(folder)))]
  d
}), fill = TRUE)
points <- read_tables("^predictions-fold[1-3]\\.csv$")
equations <- read_tables("^equations-fold[1-3]\\.csv$")
stopifnot(nrow(points) > 0, nrow(equations) > 0)
metrics <- points[, .(rmse = sqrt(sum(probability * (estimate - true)^2) / sum(probability)),
  signed_error = sum(probability * (estimate - true)) / sum(probability),
  max_absolute_error = max(abs(estimate - true)), combinations = .N,
  unique_solution = all(unique_solution),
  truncation_probability = if (kind[1L] == "beta")
    sum(probability[estimate != applied_estimate]) / sum(probability) else 0,
  max_truncation_change = if (kind[1L] == "beta") max(abs(estimate - applied_estimate)) else 0),
  by = .(mechanism, n, kind, horizon, current_R, fold, candidate)]
equation_metrics <- equations[, .(equation_rmse = sqrt(sum(probability * (estimate - true)^2) / sum(probability)),
  equation_max_error = max(abs(estimate - true))),
  by = .(mechanism, n, kind, horizon, current_R, fold, candidate)]
metrics <- merge(metrics, equation_metrics,
  by = c("mechanism", "n", "kind", "horizon", "current_R", "fold", "candidate"))
fwrite(metrics, file.path(root, "function-and-equation-errors.csv"))
fwrite(points, file.path(root, "all-function-predictions.csv"))
fwrite(equations, file.path(root, "all-equation-predictions.csv"))
theme <- theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
  plot.caption = element_text(hjust = 0, size = 11), legend.position = "none")
make_plot <- function(d, title, subtitle, xlabel, ylabel, filename, central = FALSE) {
  d <- copy(d)
  d[, sample := factor(paste0("n = ", format(n, big.mark = ",", scientific = FALSE, trim = TRUE)),
    levels = c("n = 4,000", "n = 20,000"))]
  d[, training_fit := paste("Training fit", fold)]
  limits <- range(c(d$true, d$estimate), finite = TRUE)
  padding <- max(diff(limits) * .03, .05)
  limits <- limits + c(-padding, padding)
  caption <- "Blue line: y = x. Each point is one reachable predictor combination; all points have equal size."
  if (central) {
    limits <- c(-1, 8)
    mass <- d[, .(outside = sum(probability[estimate < -1 | estimate > 8]) / sum(probability)),
      by = .(sample, training_fit)]
    caption <- paste0("View restricted to [-1, 8]. See the companion full-range plot for every point.\n",
      "Largest probability outside this view in a panel: ", sprintf("%.2f%%", 100 * max(mass$outside)), ". Blue line: y = x.")
  }
  plot <- ggplot(d, aes(true, estimate)) + geom_abline(slope = 1, intercept = 0,
    color = "#316bb0", linewidth = .8) + geom_point(color = "#24504b", size = 1.8, alpha = .72) +
    facet_grid(sample ~ training_fit) + coord_equal(xlim = limits, ylim = limits) +
    labs(title = title, subtitle = subtitle, x = xlabel, y = ylabel, caption = caption) + theme
  ggsave(file.path(root, paste0(filename, ".png")), plot, width = 12, height = 8.5, dpi = 180)
  ggsave(file.path(root, paste0(filename, ".pdf")), plot, width = 12, height = 8.5, device = grDevices::pdf)
}
for (m in c("binary_longitudinal", "discrete_dose")) {
  label <- if (m == "binary_longitudinal") "Binary treatment" else "Numerical dose"
  for (candidate_name in c("ensemble", "sieve_md")) {
    d <- points[mechanism == m & kind == "beta" & horizon == 2 & current_R == 1 & candidate == candidate_name]
    if (nrow(d)) {
      stem <- paste0(m, "-bridge-", candidate_name)
      make_plot(d, paste0(label, ": ", candidate_name, " bridge"),
        "Outcome at time 3; R_2 = 1. Raw predictions from the fitted parameterization. Same replication seed and saved splits.",
        expression("True " * beta[2](H[2], C[3])), expression("Estimated " * beta[2](H[2], C[3])), stem)
      make_plot(d, paste0(label, ": ", candidate_name, " bridge (restricted view)"),
        "The companion full-range plot retains points outside this window.",
        expression("True " * beta[2](H[2], C[3])), expression("Estimated " * beta[2](H[2], C[3])), paste0(stem, "-central"), TRUE)
    }
    if (m == "binary_longitudinal") {
      d <- points[mechanism == m & kind == "adjoint" & horizon == 2 & current_R == 1 & candidate == candidate_name]
      if (nrow(d)) make_plot(d, paste0(label, ": ", candidate_name, " adjoint"),
        "Outcome at time 3; R_2 = 1. The adjoint is unique on the reachable support. Estimated out-of-sample ratios.",
        expression("True " * lambda[2](A[2], H[2])), expression("Estimated " * lambda[2](A[2], H[2])),
        paste0(m, "-adjoint-", candidate_name))
    }
    d <- equations[mechanism == m & kind == "adjoint" & horizon == 2 & current_R == 1 & candidate == candidate_name]
    if (nrow(d)) make_plot(d, paste0(label, ": ", candidate_name, " adjoint equation"),
      "Outcome at time 3; R_2 = 1. Conditional mean of the fitted adjoint versus the true equation target.",
      "True right-hand side of the adjoint equation", "Conditional mean of fitted adjoint",
      paste0(m, "-adjoint-equation-", candidate_name))
  }
}
cat("Saved function/equation plots and metrics. Rows:", nrow(points), "\n")
