# A bridge-only wiring check on the saved n=4000 replication, not a study run.
source(file.path(Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected"), "study.R"))
name <- commandArgs(TRUE)[1L]
root <- Sys.getenv("SIM_OUTPUT")
reference <- Sys.getenv("SIM_REFERENCE")
seed <- if (name == "binary_longitudinal") 5103006L else 5203006L
saved <- readRDS(file.path(reference, paste0(name, "-n4000-seed", seed), "result.rds"))
g <- make_mechanism(name, dose_max = 3L, visits = "missing_both")
task <- study_task(saved$data, g, saved$folds, saved$learner_folds, learner_groups = 3L)$task
rows <- task$folds[[1]]$training_set; valid <- task$folds[[1]]$validation_set
args <- study_arguments(g); H <- task$vars$history("A", 2L); A <- task$vars$A[[2L]]
B <- as.matrix(task$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(task$natural[, H, drop = FALSE], saved$data[, args$health[[2L]], drop = FALSE]))
M <- saved$data[[args$measurement[2L]]]
split <- lmtp:::make_bridge_nested_folds(task$id[rows], saved$data[rows, args$measurement], task$learner_folds[[1]])
lib <- cm_library(bridge = TRUE)
lib$sieve_md$control$penalty_ids <- task$id[rows]
lib$sieve_md$control$penalty_folds <- split
started <- Sys.time()
fit <- cmbridge::fit_bridge_ensemble(B[rows, ], V[rows, ], M[rows], lib,
  fold_id = task$learner_folds[[1]], kernel = "cell")
candidate <- fit$candidates$sieve_md
stopifnot(candidate$tuning$lambda == candidate$penalty_cv$scale[candidate$penalty_cv$selected],
  candidate$tuning$lambda > 0, identical(candidate$penalty_ids, task$id[rows]),
  identical(candidate$penalty_fold_id, split(task$id[rows])),
  length(intersect(candidate$penalty_ids, task$id[valid])) == 0,
  fit$candidates$landweber$target_spec$type == "quadratic")
index <- which(candidate$penalty_cv$selected)
stopifnot(index > 1, index < nrow(candidate$penalty_cv))
grids <- list(cbind(as.data.table(candidate$penalty_cv), data.table(learner_fold = 0L)))
for (j in seq_along(fit$fold_penalty_cv)) {
  cv <- fit$fold_penalty_cv[[j]]$sieve_md
  selected <- which(cv$selected)
  stopifnot(length(selected) == 1, selected > 1, selected < nrow(cv))
  grids[[j + 1L]] <- cbind(as.data.table(cv), data.table(learner_fold = j))
}
out <- file.path(root, name); dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(grids), file.path(out, "sieve-penalty-cv.csv"))
saveRDS(fit, file.path(out, "bridge-fit.rds"))
writeLines(c(paste("Selected final ridge lambda:", candidate$tuning$lambda),
  paste("CV status:", candidate$penalty_status),
  paste("Seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))),
  "Shared splits verified; no outer validation IDs used for penalty selection."), file.path(out, "STATUS.txt"))
cat(readLines(file.path(out, "STATUS.txt")), sep = "\n")
