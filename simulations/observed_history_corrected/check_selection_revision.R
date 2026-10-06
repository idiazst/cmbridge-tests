# Four isolated checks requested after the scorer/penalty/library revision.
# The repeated-sample simulation remains stopped.
source(file.path(Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected"), "study.R"))
args <- commandArgs(TRUE)
mechanism <- args[1L]; sample_n <- as.integer(args[2L])
seed <- if (mechanism == "binary_longitudinal") 5103006L else 5203006L
root <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/selection-revision")
out <- file.path(root, paste0(mechanism, "-n", sample_n, "-seed", seed))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(warn = 1)
started <- Sys.time()
reference <- Sys.getenv("SIM_REFERENCE")
old <- if (nzchar(reference)) readRDS(file.path(reference,
  paste0(mechanism, "-n", sample_n, "-seed", seed), "result.rds")) else if (sample_n == 4000L) readRDS(file.path("results/observed-history-corrected/missing-both-study-v2/jobs",
  paste0(mechanism, "-n4000-r006.rds"))) else NULL
capture <- function(fitted, task, p, pt, g, fold) {
  saveRDS(fitted, file.path(out, paste0("nuisance-fits-fold", fold, ".rds")))
  true <- population_truth_functions(p, g)
  points <- equations <- selection <- boundaries <- scores <- list()
  for (s in seq_len(g$tau)) {
    H <- pt$vars$history("A", s); A <- pt$vars$A[[s]]
    B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
    V <- as.matrix(cbind(pt$natural[, H, drop = FALSE], p[, study_arguments(g)$health[[s]], drop = FALSE]))
    observed <- p[[paste0("R", s + 1L)]] == 1
    current <- p[[paste0("R", s)]]
    bkey <- cell_key(B); vkey <- rep(NA_character_, nrow(p))
    vkey[observed] <- cell_key(V[observed, , drop = FALSE])
    rhs <- conditional_mean(apply(true$ratios[observed, seq_len(s), drop = FALSE], 1, prod) *
      p[[paste0("Y", s + 1L)]][observed], w = p$probability[observed], key = vkey[observed])
    for (kind in c("beta", "adjoint")) {
      ensemble <- fitted[[paste0(kind, "_fits")]][[s]]
      models <- c(list(ensemble = ensemble), ensemble$candidates)
      scoring <- attr(ensemble$fold_gram[[1L]], "cell_scoring")
      scores[[length(scores) + 1L]] <- data.table(kind = kind, horizon = s, fold = fold,
        candidate = names(ensemble$weights), weight = as.numeric(ensemble$weights),
        raw_u_loss = as.numeric(ensemble$candidate_raw_cv_loss),
        projected_loss = as.numeric(ensemble$candidate_cv_loss),
        projection_norm = ensemble$psd_projection$adjustment_norm,
        negative_eigenvalues = ensemble$psd_projection$removed_negative_eigenvalues)
      for (j in seq_along(ensemble$fold_gram)) {
        metadata <- attr(ensemble$fold_gram[[j]], "cell_scoring")
        for (name in names(ensemble$fold_penalty_cv[[j]])) {
          cv <- ensemble$fold_penalty_cv[[j]][[name]]
          if (!is.null(cv)) selection[[length(selection) + 1L]] <- cbind(as.data.table(cv),
            data.table(kind = kind, horizon = s, fold = fold, learner_fold = j, candidate = name))
          trace <- ensemble$fold_penalty_boundary[[j]][[name]]
          if (!is.null(trace) && nrow(trace)) boundaries[[length(boundaries) + 1L]] <- cbind(as.data.table(trace),
            data.table(kind = kind, horizon = s, fold = fold, learner_fold = j, candidate = name))
        }
      }
      for (name in names(models)) {
        model <- models[[name]]
        estimate <- rep(0, nrow(p))
        if (kind == "beta") estimate[observed] <- predict(model, V[observed, , drop = FALSE]) else
          estimate <- predict(model, B)
        rows <- if (kind == "beta") which(observed) else seq_len(nrow(p))
        table <- data.table(cell = if (kind == "beta") vkey[rows] else bkey[rows],
          current_R = current[rows], true = true[[if (kind == "beta") "beta" else "lambda"]][rows, s],
          estimate = estimate[rows], probability = p$probability[rows])
        table <- table[, .(current_R = current_R[1L], true = true[1L], estimate = estimate[1L],
                          probability = sum(probability)), by = cell]
        table[, `:=`(kind = kind, horizon = s, fold = fold, candidate = name,
          unique_solution = if (kind == "beta") current_R == 1 else !g$dose | current_R == 0)]
        if (kind == "beta") table[, applied_estimate := pmin(100, pmax(-100, estimate))]
        points[[length(points) + 1L]] <- table
        if (kind == "beta") {
          lhs <- conditional_mean(estimate, w = p$probability, key = bkey)
          q <- data.table(cell = names(lhs), estimate = as.numeric(lhs), true = 1)
          reference <- data.table(cell = bkey, probability = p$probability, current_R = current)
        } else {
          lhs <- conditional_mean(estimate[observed], w = p$probability[observed], key = vkey[observed])
          q <- data.table(cell = names(lhs), estimate = as.numeric(lhs), true = as.numeric(rhs[names(lhs)]))
          reference <- data.table(cell = vkey[observed], probability = p$probability[observed], current_R = current[observed])
        }
        reference <- reference[, .(probability = sum(probability), current_R = current_R[1L]), by = cell]
        q <- merge(q, reference, by = "cell")
        q[, `:=`(kind = kind, horizon = s, fold = fold, candidate = name)]
        equations[[length(equations) + 1L]] <- q
        if (!is.null(model$penalty_cv)) selection[[length(selection) + 1L]] <- cbind(
          as.data.table(model$penalty_cv), data.table(kind = kind, horizon = s, fold = fold,
            learner_fold = 0L, candidate = name))
        if (!is.null(model$penalty_boundary) && nrow(model$penalty_boundary)) boundaries[[length(boundaries) + 1L]] <- cbind(
          as.data.table(model$penalty_boundary), data.table(kind = kind, horizon = s, fold = fold,
            learner_fold = 0L, candidate = name))
      }
    }
  }
  for (pair in list(list(points, "predictions"), list(equations, "equations"),
                    list(selection, "penalties"), list(boundaries, "boundary-attempts"), list(scores, "weights")))
    fwrite(rbindlist(pair[[1L]], fill = TRUE), file.path(out, paste0(pair[[2L]], "-fold", fold, ".csv")))
  cat(format(Sys.time()), "Saved nuisance diagnostics for outer split", fold, "\n")
}
if (Sys.getenv("SIM_EXPORT_ONLY") == "1") {
  set.seed(seed); g <- make_mechanism(mechanism, dose_max = 3L, visits = "missing_both")
  d <- draw_data(sample_n, g)
  task <- study_task(d, g, folds = if (is.null(old)) 3L else old$folds,
    learner_folds = if (is.null(old)) NULL else old$learner_folds,
    learner_groups = 3L, balance_measurement = TRUE)$task
  p <- enumerate_data(g); pt <- study_task(p, g)$task
  for (fold in 1:3) {
    file <- file.path(out, paste0("nuisance-fits-fold", fold, ".rds"))
    if (file.exists(file)) capture(readRDS(file), task, p, pt, g, fold)
  }
  quit(save = "no")
}
cache <- function(fold, constant) {
  file <- file.path(out, paste0("nuisance-fits-fold", fold, ".rds"))
  if (!constant && file.exists(file)) return(readRDS(file))
  NULL
}
tryCatch({
  result <- run_dataset(mechanism, sample_n, 6L, seed, all_correct = TRUE,
    dose_max = 3L, learner_groups = 3L,
    balance_measurement = TRUE, folds = if (is.null(old)) 3L else old$folds,
    learner_folds = if (is.null(old)) NULL else old$learner_folds,
    visits = "missing_both", store_scores = TRUE, nuisance_callback = capture, nuisance_cache = cache)
  if (!is.null(old)) {
    stopifnot(identical(result$folds, old$folds), identical(result$learner_folds, old$learner_folds))
    if (!is.null(old$data)) stopifnot(identical(result$data, old$data))
  }
  saveRDS(result, file.path(out, "result.rds"))
  fwrite(result$results, file.path(out, "estimates.csv"))
  fwrite(result$errors, file.path(out, "estimator-errors.csv"))
  fwrite(result$tuning, file.path(out, "estimator-penalties.csv"))
  fwrite(result$ensemble_weights, file.path(out, "all-weights.csv"))
  fwrite(result$eif_components, file.path(out, "eif-components.csv"))
  fwrite(result$candidate_failures, file.path(out, "nested-candidate-failures.csv"))
  fwrite(result$penalty_trial_failures, file.path(out, "nested-penalty-trial-failures.csv"))
  writeLines(c(paste("Completed", Sys.time()), paste("Seconds", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    paste("Estimator error rows", nrow(result$errors))), file.path(out, "STATUS.txt"))
}, error = function(e) {
  saveRDS(e, file.path(out, "failure.rds"))
  if (!is.null(e$penalty_cv)) fwrite(e$penalty_cv, file.path(out, "failed-penalty-grid.csv"))
  if (!is.null(e$boundary_history)) fwrite(e$boundary_history, file.path(out, "failed-boundary-attempts.csv"))
  writeLines(c(paste("Failed", Sys.time()), conditionMessage(e)), file.path(out, "STATUS.txt"))
  stop(e)
})
