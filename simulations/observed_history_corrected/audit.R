source("simulations/observed_history_corrected/dgp.R")
out <- "results/observed-history-corrected"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- list()
truth <- list()
for (name in c("binary_point", "binary_longitudinal", "discrete_dose")) {
  mechanism <- make_mechanism(name, visits = "scheduled")
  audit <- require_valid_mechanism(mechanism)
  stopifnot(max(abs(exact_truth(mechanism) - backward_truth(mechanism))) < 1e-12)
  rows[[name]] <- audit
  truth[[name]] <- data.frame(mechanism = name, outcome_time = seq_len(mechanism$tau) + 1L,
                             truth = exact_truth(mechanism),
                             backward_truth = backward_truth(mechanism))
  set.seed(20261004)
  pilot <- draw_data(1000, mechanism)
  if (mechanism$tau == 2L) stopifnot(all(pilot$R2 == 1))
  stopifnot(all(is.finite(unlist(pilot[paste0("R", seq_len(mechanism$tau + 1L))]))))
}
write.csv(do.call(rbind, rows), file.path(out, "assumption_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, truth), file.path(out, "truth.csv"), row.names = FALSE)
# Demonstrate why imposing deterministic continuation alone is insufficient.
failed <- make_mechanism("binary_longitudinal", visits = "intermittent")
audit <- audit_mechanism(failed)
stopifnot(all(audit$assumption_1), all(audit$assumption_3), all(audit$assumption_5),
          any(!audit$assumption_7))
rejection <- tryCatch({draw_data(100, failed); FALSE}, error = function(e) {
  grepl("assumption_7", conditionMessage(e), fixed = TRUE)
})
stopifnot(rejection)
write.csv(audit, file.path(out, "rejected_intermittent_design.csv"), row.names = FALSE)
cat("All three scheduled-treatment designs pass the finite-support checks.\n")
cat("Independent treatment/health/measurement random draws establish Assumption 2;\n")
cat("measurement uses only H_t and C_(t+1), establishing Assumption 4 by construction.\n")
cat("The intermittent design with deterministic treatment passes 1/3/5 but fails 7 and is blocked.\n")
cat("Only audit outputs and 1,000-person generation pilots were produced; no estimator study has run.\n")
