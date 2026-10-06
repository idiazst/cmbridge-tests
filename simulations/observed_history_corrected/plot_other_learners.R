# Compare saved candidates without rerunning any nuisance fits.
library(data.table)
library(ggplot2)
root <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/selection-revision")
points <- fread(file.path(root, "all-function-predictions.csv"))
equations <- fread(file.path(root, "all-equation-predictions.csv"))
candidate_labels <- c(landweber = "Landweber", pmmr = "PMMR",
  saturated_cv = "Saturated L1", ensemble = "Ensemble")
make_plot <- function(d, title, subtitle, xlabel, ylabel, stem) {
  d <- copy(d)
  d[, learner := factor(candidate, levels = names(candidate_labels), labels = candidate_labels)]
  d[, training_fit := paste("Training fit", fold)]
  limits <- range(c(d$true, d$estimate), finite = TRUE)
  padding <- max(diff(limits) * .03, .05)
  limits <- limits + c(-padding, padding)
  p <- ggplot(d, aes(true, estimate)) +
    geom_abline(slope = 1, intercept = 0, color = "#316bb0", linewidth = .7) +
    geom_point(color = "#24504b", alpha = .7, size = 1.4) +
    facet_grid(learner ~ training_fit, drop = TRUE) +
    coord_equal(xlim = limits, ylim = limits) +
    labs(title = title, subtitle = subtitle, x = xlabel, y = ylabel,
      caption = "Blue line: y = x. All reachable combinations retained, with equal point sizes. Saved fits; no refitting.") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
  height <- 1.5 + 2.6 * length(unique(d$candidate))
  ggsave(file.path(root, paste0(stem, ".png")), p, width = 12, height = height, dpi = 160)
  ggsave(file.path(root, paste0(stem, ".pdf")), p, width = 12, height = height, device = grDevices::pdf)
}
for (mechanism_name in c("binary_longitudinal", "discrete_dose")) {
  label <- if (mechanism_name == "binary_longitudinal") "Binary treatment" else "Numerical dose"
  for (sample_n in sort(unique(points$n))) {
    subtitle <- paste0("n = ", format(sample_n, big.mark = ",", trim = TRUE),
      "; outcome at time 3; R_2 = 1. Raw fitted predictions.")
    selected <- points[mechanism == mechanism_name & n == sample_n &
      kind == "beta" & horizon == 2 & current_R == 1 & candidate %in% names(candidate_labels)]
    make_plot(selected, paste0(label, ": other bridge learners and ensemble"), subtitle,
      expression("True " * beta[2](H[2], C[3])), expression("Estimated " * beta[2](H[2], C[3])),
      paste0(mechanism_name, "-n", sample_n, "-other-bridge-learners"))
    selected <- equations[mechanism == mechanism_name & n == sample_n &
      kind == "adjoint" & horizon == 2 & current_R == 1 & candidate %in% names(candidate_labels)]
    make_plot(selected, paste0(label, ": other adjoint learners and ensemble"),
      sub("Raw fitted predictions.", "Defining conditional equation.", subtitle, fixed = TRUE),
      "True right-hand side of the adjoint equation", "Conditional mean of fitted adjoint",
      paste0(mechanism_name, "-n", sample_n, "-other-adjoint-learners"))
  }
}
cat("Saved other-learner comparisons for all four datasets.\n")
