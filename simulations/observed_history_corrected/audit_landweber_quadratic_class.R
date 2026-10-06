# Linear-algebra membership audit only. No fitted oracle functions or new sample.
source(file.path(Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected"), "study.R"))
output <- Sys.getenv("SIM_OUTPUT")
rows <- list()
for (name in c("binary_longitudinal", "discrete_dose")) {
  g <- make_mechanism(name, dose_max = 3L, visits = "missing_both")
  require_valid_mechanism(g)
  population <- enumerate_data(g); task <- study_task(population, g)$task
  truth <- population_truth_functions(population, g)
  for (s in 1:2) {
    H <- task$vars$history("A", s)
    observed <- population[[paste0("R", s + 1L)]] == 1
    x <- as.matrix(cbind(task$natural[observed, H, drop = FALSE],
      population[observed, study_arguments(g)$health[[s]], drop = FALSE]))
    key <- cell_key(x); keep <- !duplicated(key)
    beta <- truth$beta[observed, s][keep]
    current_R <- population[[paste0("R", s)]][observed][keep]
    x <- x[keep, , drop = FALSE]
    control <- cm_library(bridge = TRUE)$landweber$control
    stopifnot(control$target_basis == "quadratic", control$link == "inverse_logit")
    spec <- cmbridge:::.fit_basis_spec(x, control$target_basis)
    stopifnot(identical(spec, cmbridge:::.fit_basis_spec(x[1, , drop = FALSE], control$target_basis)))
    design <- cmbridge:::.eval_basis_spec(x, spec)
    eta <- -log(beta - 1)
    decomposition <- svd(design)
    retained <- decomposition$d > max(decomposition$d) * 1e-10
    coefficients <- as.numeric(decomposition$v[, retained, drop = FALSE] %*%
      (as.numeric(crossprod(decomposition$u[, retained, drop = FALSE], eta)) / decomposition$d[retained]))
    reproduced <- 1 + exp(-as.numeric(design %*% coefficients))
    stopifnot(all(is.finite(coefficients)), max(abs(reproduced - beta)) < 1e-10)
    for (state in unique(current_R)) {
      ii <- current_R == state
      rows[[length(rows) + 1L]] <- data.table(mechanism = name, horizon = s,
        current_R = state, reachable_combinations = sum(ii), basis_columns = ncol(design),
        maximum_bridge_representation_error = max(abs(reproduced[ii] - beta[ii])),
        maximum_link_representation_error = max(abs(as.numeric(design %*% coefficients)[ii] - eta[ii])),
        basis_independent_of_training_combinations = TRUE)
    }
  }
}
audit <- rbindlist(rows)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(audit, file.path(output, "landweber-function-class.csv"))
print(audit)
cat("Verified true bridge membership on the entire reachable support, including missed visits.\n")
