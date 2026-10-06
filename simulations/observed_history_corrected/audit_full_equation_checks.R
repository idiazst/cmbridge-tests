# Separate read-only audit for full-equation cell-scored checks.
source(file.path(Sys.getenv("STUDY_SOURCE"), "study.R"))
root <- Sys.getenv("SIM_OUTPUT"); reference <- Sys.getenv("SIM_REFERENCE")
rows <- list(); solver_rows <- list()
for (name in c("binary_longitudinal", "discrete_dose")) {
  seed <- if (name == "binary_longitudinal") 5103006L else 5203006L
  folder <- file.path(root, paste0(name, "-n4000-seed", seed))
  result <- readRDS(file.path(folder, "result.rds"))
  old <- readRDS(file.path(reference, basename(folder), "result.rds"))
  stopifnot(nrow(result$errors) == 0, nrow(result$results) == 4,
    identical(result$data, old$data), identical(result$folds, old$folds),
    identical(result$learner_folds, old$learner_folds))
  g <- make_mechanism(name, dose_max = 3L, visits = "missing_both")
  task <- study_task(result$data, g, result$folds, result$learner_folds, learner_groups = 3L)$task
  args <- study_arguments(g)
  for (fold in 1:3) {
    fits <- readRDS(file.path(folder, paste0("nuisance-fits-fold", fold, ".rds")))
    training <- task$folds[[fold]]$training_set; validation <- task$folds[[fold]]$validation_set
    split <- lmtp:::make_bridge_nested_folds(task$id[training], result$data[training, args$measurement],
      task$learner_folds[[fold]])
    for (s in 1:2) for (kind in c("beta", "adjoint")) {
      model <- fits[[paste0(kind, "_fits")]][[s]]
      stopifnot(identical(names(model$weights), c("sieve_md", "landweber", "pmmr")))
      for (candidate in model$candidates) {
        stopifnot(candidate$tuning$link == if (kind == "beta") "inverse_logit" else "identity",
          all(is.finite(fitted(candidate)[!is.na(fitted(candidate))])))
        if (kind == "beta") stopifnot(all(fitted(candidate)[!is.na(fitted(candidate))] >= 1))
      }
      if (kind == "beta") stopifnot(model$candidates$landweber$target_spec$type == "poly",
        model$candidates$landweber$target_spec$degree == 1L,
        model$candidates$landweber$instrument_spec$type == "cell",
        model$candidates$sieve_md$instrument_spec$type == "cell",
        model$candidates$sieve_md$target_spec$type == "cell_linear")
      for (candidate in model$candidates) {
        if (!is.null(candidate$solver)) solver_rows[[length(solver_rows)+1L]] <- data.table(
          mechanism=name,fold=fold,horizon=s,kind=kind,candidate=candidate$method,
          tolerance_reached=if(candidate$method=="landweber")isTRUE(candidate$solver$stopped_by_tolerance)
            else isTRUE(candidate$solver$converged),
          gradient=candidate$solver$coefficient_gradient)
      }
      stopifnot(length(model$candidate_failures) == 0L)
      sieve <- model$candidates$sieve_md
      for (candidate in model$candidates) if(!is.null(candidate$penalty_cv)) {
        cv_candidate <- candidate$penalty_cv; chosen <- which(cv_candidate$selected)
        stopifnot(length(chosen)==1L,chosen>1L,chosen<nrow(cv_candidate),
          all(cv_candidate$scale>0),!cv_candidate$failed[chosen],
          is.finite(cv_candidate$loss[chosen]),cv_candidate$loss[chosen]==min(cv_candidate$loss),
          candidate$tuning[[candidate$penalty_parameter]]==cv_candidate$scale[chosen],
          identical(candidate$penalty_ids,task$id[training]),
          identical(candidate$penalty_fold_id,split(task$id[training])),
          length(intersect(candidate$penalty_ids,task$id[validation]))==0L)
      }
      cv <- sieve$penalty_cv; selected <- which(cv$selected)
      stopifnot(length(selected) == 1, selected > 1, selected < nrow(cv), all(cv$scale > 0),
        sieve$tuning$lambda == cv$scale[selected],
        identical(sieve$penalty_ids, task$id[training]),
        identical(sieve$penalty_fold_id, split(task$id[training])),
        length(intersect(sieve$penalty_ids, task$id[validation])) == 0)
      rows[[length(rows) + 1L]] <- data.table(mechanism = name, n = 4000L, fold, horizon = s,
        kind, lambda = sieve$tuning$lambda, status = sieve$penalty_status,
        grid_values = nrow(cv), boundary_failures = nrow(sieve$penalty_boundary),
        minimum_raw_u = cv$loss[selected], penalty_training_n = length(sieve$penalty_ids))
      for (cv_fold in model$fold_penalty_cv) {
        for(path in cv_fold) if(!is.null(path)) {
          selected <- which(path$selected)
          stopifnot(length(selected)==1L,selected>1L,selected<nrow(path),all(path$scale>0),
            !path$failed[selected],is.finite(path$loss[selected]),path$loss[selected]==min(path$loss))
        }
      }
    }
  }
  stopifnot(max(abs(result$results$mean_sequential_eif + result$results$mean_bridge_adjoint_eif -
    (result$results$estimate - result$results$truth))) < 1e-10)
}
fwrite(rbindlist(rows), file.path(root, "sieve-selected-penalties.csv"))
solver_table<-rbindlist(solver_rows)
fwrite(solver_table,file.path(root,"solver-tolerance-audit.csv"))
stopifnot(all(solver_table$tolerance_reached))
cat("Verified CV-selected bridge/adjoint ridge penalties, positive interior selections, shared training-only splits, unchanged data and outer folds, all bridge links, and EIF sums.\n")
