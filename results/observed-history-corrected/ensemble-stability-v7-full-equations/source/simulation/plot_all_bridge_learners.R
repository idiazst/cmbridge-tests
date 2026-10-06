library(data.table)
root <- Sys.getenv("SIM_OUTPUT")
points <- fread(file.path(root, "all-function-predictions.csv"))
labels <- c(sieve_md = "Sieve", landweber = "Landweber", pmmr = "PMMR", ensemble = "Ensemble")
colors <- c("#2166ac", "#b35806", "#5e3c99")
for (mechanism_name in c("binary_longitudinal", "discrete_dose")) {
  d <- points[mechanism == mechanism_name & kind == "beta" & horizon == 2 & current_R == 1]
  stopifnot(nrow(d) > 0, all(d$unique_solution), all(is.finite(d$estimate)), all(d$estimate >= 1))
  title <- if (mechanism_name == "binary_longitudinal") "Binary treatment" else "Numerical dose"
  sample_sizes <- sort(unique(d$n))
  draw <- function() {
    par(pty = "s", mfrow = if (length(sample_sizes) == 1L) c(2, 2) else c(4, length(sample_sizes)),
      mar = c(3.7, 4.0, 2.3, .8), oma = c(3.5, .5, 3.0, .5))
    for (candidate_name in names(labels)) {
      limits <- range(d[candidate == candidate_name, c(true, estimate)])
      padding <- max(diff(limits) * .05, .1); limits <- limits + c(-padding, padding)
      for (sample_n in sample_sizes) {
        z <- d[candidate == candidate_name & n == sample_n]
        plot(z$true, z$estimate, xlim = limits, ylim = limits, asp = 1, pch = 19, cex = .55,
          col = adjustcolor(colors[z$fold], alpha.f = .55),
          xlab = expression("True " * beta[2]), ylab = expression("Estimated " * beta[2]),
          main = "")
        graphics::title(main = paste0(labels[candidate_name], "; n = ", format(sample_n, big.mark = ",")), line = .5, cex.main = .95)
        abline(0, 1, lty = 2, lwd = 1.5, col = "#333333")
      }
    }
    mtext(paste0(title, ": all bridge candidates use inverse-expit"), outer = TRUE, side = 3, line = 1, cex = 1.1, font = 2)
    mtext("Outcome at time 3; R_2 = 1. All reachable combinations and all three training fits shown.",
      outer = TRUE, side = 1, line = .7, cex = .78)
    mtext("Blue: fit 1. Orange: fit 2. Purple: fit 3. Dashed line: y = x. Equal axes within each panel.",
      outer = TRUE, side = 1, line = 1.8, cex = .78)
  }
  stem <- file.path(root, paste0(mechanism_name, "-all-bridge-learners"))
  height <- if (length(sample_sizes) == 1L) 10 else 13
  pdf(paste0(stem, ".pdf"), width = 10, height = height); draw(); dev.off()
  png(paste0(stem, ".png"), width = 1800, height = height*180, res = 180); draw(); dev.off()
}
cat("Saved complete bridge-candidate comparisons for the available sample sizes.\n")
