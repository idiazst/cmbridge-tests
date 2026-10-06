library(data.table)
library(cmbridge)
root <- Sys.getenv("SIM_OUTPUT"); reference <- Sys.getenv("SIM_REFERENCE")
rows <- list()
for (mechanism in c("binary_longitudinal", "discrete_dose")) for (sample_n in c(4000L, 20000L)) {
  seed <- if (mechanism == "binary_longitudinal") 5103006L else 5203006L
  folder <- file.path(root, paste0(mechanism, "-n", sample_n, "-seed", seed))
  result <- readRDS(file.path(folder, "result.rds"))
  old_result <- readRDS(file.path(reference, basename(folder), "result.rds"))
  stopifnot(nrow(result$errors) == 0, nrow(result$results) == 4,
    identical(result$data, old_result$data), identical(result$folds, old_result$folds),
    identical(result$learner_folds, old_result$learner_folds))
  for (fold in 1:3) {
    fits <- readRDS(file.path(folder, paste0("nuisance-fits-fold", fold, ".rds")))
    for (s in 1:2) for (kind in c("beta", "adjoint")) {
      model <- fits[[paste0(kind, "_fits")]][[s]]
      stopifnot(identical(names(model$weights), c("sieve_md", "landweber", "pmmr")),
        all(model$weights >= 0), abs(sum(model$weights) - 1) < 1e-10)
      for (name in names(model$candidates)) {
        fit <- model$candidates[[name]]
        values <- fitted(fit); values <- values[!is.na(values)]
        expected_link <- if (kind == "beta") "inverse_logit" else "identity"
        stopifnot(fit$tuning$link == expected_link, all(is.finite(values)))
        penalty <- NA_real_
        if (kind == "beta") {
          stopifnot(all(values >= 1), all(is.finite(fit$link_coefficients)))
          if (name == "sieve_md") {
            penalty <- fit$tuning$lambda * sum(fit$coefficients[-1L]^2)
            stopifnot(fit$solver$converged,
              max(fit$solver$coefficient_gradient, fit$solver$function_gradient) <= fit$solver$tolerance,
              abs(fit$solver$penalized_loss - (fit$moment_loss + penalty)) < 1e-7)
          } else if (name == "landweber") {
            penalty <- 0
            stopifnot(fit$solver$final_loss <= fit$solver$initial_loss + 1e-12,
              abs(fit$solver$final_loss - fit$moment_loss) < 1e-10)
          } else {
            penalty <- fit$tuning$lambda * sum(fit$coefficients^2)
            stopifnot(fit$solver$converged,
              fit$solver$coefficient_gradient <= fit$solver$tolerance,
              abs(fit$solver$penalized_loss - (fit$moment_loss + penalty)) < 1e-10)
          }
        }
        rows[[length(rows) + 1L]] <- data.table(mechanism, n = sample_n, fold, horizon = s,
          kind, candidate = name, link = expected_link, minimum_value = min(values), maximum_value = max(values),
          coefficient_gradient = if (kind == "beta") fit$solver$coefficient_gradient else NA_real_,
          moment_loss = fit$moment_loss, penalty)
      }
    }
  }
}
fwrite(rbindlist(rows), file.path(root, "all-bridge-link-audit.csv"))
new <- fread(file.path(root, "all-function-predictions.csv"))
old <- fread(file.path(reference, "all-function-predictions.csv"))
keys <- c("mechanism", "n", "kind", "horizon", "current_R", "fold", "candidate", "cell")
a <- new[kind == "adjoint"]; b <- old[kind == "adjoint"]
comparison <- merge(a[, c(keys, "estimate"), with = FALSE], b[, c(keys, "estimate"), with = FALSE],
  by = keys, suffixes = c("_new", "_old"))
stopifnot(nrow(comparison) == nrow(a), max(abs(comparison$estimate_new - comparison$estimate_old)) < 1e-10,
  all(is.finite(new[kind == "beta"]$estimate)), all(new[kind == "beta"]$estimate >= 1))
writeLines(c(paste("Unchanged adjoint prediction rows:", nrow(comparison)),
  paste("Maximum adjoint prediction difference:", max(abs(comparison$estimate_new - comparison$estimate_old))),
  paste("Minimum raw bridge prediction:", min(new[kind == "beta"]$estimate)),
  "All four data sets, outer splits and learner-validation splits match the reference."),
  file.path(root, "all-bridge-link-audit.txt"))
cat("Verified all bridge links, objective consistency, finite predictions, unchanged adjoints, data and splits.\n")
