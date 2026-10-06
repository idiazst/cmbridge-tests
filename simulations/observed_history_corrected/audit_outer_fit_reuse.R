source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
root <- Sys.getenv('SIM_OUTPUT'); reference <- Sys.getenv('SIM_REFERENCE')
checks <- list()
for (mechanism in c('binary_longitudinal', 'discrete_dose')) {
  seed <- if (mechanism == 'binary_longitudinal') 5103006L else 5203006L
  folder <- paste0(mechanism, '-n4000-seed', seed)
  original <- readRDS(file.path(reference, folder, 'result.rds'))
  set.seed(seed); g <- make_mechanism(mechanism, dose_max = 3L, visits = 'missing_both')
  d <- draw_data(4000L, g)
  prepared <- study_task(d, g, original$folds, original$learner_folds, learner_groups = 3L)
  task <- prepared$task; args <- study_arguments(g); fold <- 2L; horizon <- 2L
  rows <- task$folds[[fold]]$training_set
  nested <- lmtp:::make_bridge_nested_folds(task$id[rows], d[rows, args$measurement, drop = FALSE],
    task$learner_folds[[fold]])
  library <- cm_library(FALSE, TRUE)
  for (i in seq_along(library)) if (!is.null(library[[i]]$control$penalty_scales)) {
    library[[i]]$control$penalty_ids <- task$id[rows]
    library[[i]]$control$penalty_folds <- nested
  }
  H <- task$vars$history('A', horizon); A <- task$vars$A[[horizon]]
  B <- as.matrix(task$natural[rows, c(H, A), drop = FALSE])
  V <- as.matrix(cbind(task$natural[rows, H, drop = FALSE], d[rows, args$health[[horizon]], drop = FALSE]))
  M <- d[[args$measurement[horizon]]][rows]
  current <- cmbridge::fit_bridge_ensemble(B, V, M, library = library,
    fold_id = task$learner_folds[[fold]], kernel = study_control()$.cm_kernel)
  cached <- readRDS(file.path(root, folder, paste0('nuisance-fits-fold', fold, '.rds')))$beta_fits[[horizon]]
  prediction_difference <- max(abs(predict(current, V[M == 1, , drop = FALSE]) -
    predict(cached, V[M == 1, , drop = FALSE])))
  weight_difference <- max(abs(current$weights - cached$weights))
  lambda_difference <- abs(current$candidates$sieve_md$tuning$lambda - cached$candidates$sieve_md$tuning$lambda)
  stopifnot(prediction_difference < 1e-10, weight_difference < 1e-10, lambda_difference == 0,
    identical(current$candidates$sieve_md$penalty_fold_id, cached$candidates$sieve_md$penalty_fold_id))
  checks[[length(checks) + 1L]] <- data.table(mechanism, fold, horizon,
    prediction_difference, weight_difference, lambda_difference)
}
fwrite(rbindlist(checks), file.path(root, 'outer-fit-reuse-audit.csv'))
cat('Fresh outer bridge refits reproduce cached predictions, weights, penalties, and shared splits.\n')
