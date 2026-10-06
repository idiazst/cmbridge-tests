source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- Sys.getenv('DEBUG_OUTPUT', 'results/observed-history-corrected/ensemble-small-sample-debug')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
capture <- new.env(parent = emptyenv()); capture$count <- 0L
assign('.small_ensemble_capture', capture, .GlobalEnv)
assign('.small_ensemble_out', normalizePath(out), .GlobalEnv)
trace('.fit_ensemble', where = asNamespace('cmbridge'), print = FALSE,
 tracer = quote({
   if (length(M) <= 50L) {
     .GlobalEnv$.small_ensemble_capture$count <- .GlobalEnv$.small_ensemble_capture$count + 1L
     saveRDS(list(B = B, V = V, M = M, phi = phi, type = type, library = library,
       n_folds = n_folds, fold_id = fold_id, seed = seed, kernel = kernel,
       kernel_control = kernel_control, solver_control = solver_control),
       file.path(.GlobalEnv$.small_ensemble_out,
         sprintf('ensemble-input-%04d.rds', .GlobalEnv$.small_ensemble_capture$count)))
   }
 }))
trace('.fit_landweber', where = asNamespace('cmbridge'), print = FALSE,
 exit = quote({
   if (length(y) <= 100L) {
     key <- digest::digest(list(returnValue()$coefficients, returnValue()$target_spec))
     assign(key, list(y = y, d = d, x = x, z = z, control = control), .GlobalEnv$.small_ensemble_capture)
   }
 }))
trace('predict.cmbridge_fit', where = asNamespace('cmbridge'), print = FALSE,
 exit = quote({
   if (object$method == 'landweber' && any(!is.finite(returnValue()))) {
     key <- digest::digest(list(object$coefficients, object$target_spec))
     saveRDS(list(object = object, newdata = newdata,
       input = if (exists(key, .GlobalEnv$.small_ensemble_capture, inherits = FALSE))
         get(key, .GlobalEnv$.small_ensemble_capture, inherits = FALSE) else NULL),
       file.path(.GlobalEnv$.small_ensemble_out, 'landweber-nonfinite.rds'))
   }
 }))
result <- run_dataset('discrete_dose', 500L, 1L, 5201001L,
 all_correct = TRUE, dose_max = 3L, learner_groups = 3L, balance_measurement = TRUE)
saveRDS(result, file.path(out, 'diagnostic-result.rds'))
print(result$errors)
cat('Captured', capture$count, 'small ensemble inputs.\n')
