library(data.table)
root <- Sys.getenv("SIM_OUTPUT")
rows <- list()
for (name in c("binary_longitudinal", "discrete_dose")) {
  seed <- if (name == "binary_longitudinal") 5103006L else 5203006L
  folder <- file.path(root, paste0(name, "-n4000-seed", seed))
  for (fold in 1:3) {
    fits <- readRDS(file.path(folder, paste0("nuisance-fits-fold", fold, ".rds")))
    for (s in 1:2) {
      fit <- fits$beta_fits[[s]]$candidates$landweber
      stopifnot(fit$target_spec$type == "quadratic", fit$solver$final_loss <= fit$solver$initial_loss)
      rows[[length(rows)+1L]] <- data.table(mechanism = name, n = 4000L, fold,
        horizon = s, iterations = fit$solver$iterations, limit = fit$tuning$max_iter,
        stopped_by_tolerance = fit$solver$stopped_by_tolerance,
        initial_loss = fit$solver$initial_loss, final_loss = fit$solver$final_loss,
        coefficient_gradient = fit$solver$coefficient_gradient,
        final_delta = fit$tuning$final_delta, tolerance = fit$tuning$tol,
        backtracking_reductions = fit$solver$backtracking_reductions)
    }
  }
}
fwrite(rbindlist(rows), file.path(root, "landweber-stopping.csv"))
print(rbindlist(rows))
