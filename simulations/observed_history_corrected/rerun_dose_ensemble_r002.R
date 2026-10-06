# Refit the original n=4,000 sample with its original estimator and learner splits.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/dose-ensemble-r002')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(out, 'discrete_dose-n4000-r002.rds')
if (file.exists(checkpoint)) stop('Completed checkpoint already exists; preserve it.')
input <- 'results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds'
reference <- readRDS(input)
stopifnot(reference$design$n == 4000L, reference$design$seed == 5203002L,
          identical(vapply(cm_library(), function(x) x$method, ''),
            c(sieve_md = 'sieve_md', landweber = 'landweber', pmmr = 'pmmr', saturated_cv = 'saturated_l1')),
          'SL.earth' %in% sl_library, !('SL.earth' %in% trt_library),
          cm_library()$sieve_md$control$target_basis == 'cell',
          cm_library()$sieve_md$control$instrument_basis == 'cell')
writeLines(c(capture.output(sessionInfo()),
  paste('original input SHA256', digest::digest(file = input, algo = 'sha256')),
  'Seed: 5203002; original n=4,000 data; original estimator and learner splits.',
  'MARS only in sequential outcome regressions; four cmbridge candidates.'),
  file.path(out, 'session.txt'))
trace('bridge_nuisance_fit', where = asNamespace('lmtp'), print = FALSE,
  tracer = quote(cat(format(Sys.time()), 'bridge/adjoint nuisance sample', fold, '\n')))
trace('bridge_response_engine', where = asNamespace('lmtp'), print = FALSE,
  tracer = quote(cat(format(Sys.time()), 'sequential response engine, sample', fold, '\n')))
started <- proc.time()[3L]; warnings <- character()
result <- withCallingHandlers(run_dataset('discrete_dose', 4000L, 2L, 5203002L,
  all_correct = TRUE, dose_max = 3L, learner_groups = 3L, balance_measurement = TRUE,
  folds = reference$folds, learner_folds = reference$learner_folds),
  warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
untrace('bridge_nuisance_fit', where = asNamespace('lmtp'))
untrace('bridge_response_engine', where = asNamespace('lmtp'))
result$seconds <- proc.time()[3L] - started
result$job_warnings <- unique(warnings)
stopifnot(identical(result$folds, reference$folds),
          identical(result$learner_folds, reference$learner_folds))
saveRDS(result, paste0(checkpoint, '.tmp')); stopifnot(file.rename(paste0(checkpoint, '.tmp'), checkpoint))
for (name in c('results', 'diagnostics', 'errors', 'tuning', 'ensemble_weights', 'occupancy'))
  fwrite(result[[name]], file.path(out, paste0(name, '.csv')))
if (nrow(result$errors)) {print(result$errors); stop('Rerun retained errors; inspect the saved checkpoint.')}
old <- reference$results[beta_correct & lambda_correct & ratio_correct & m_correct]
comparison <- merge(old[, .(horizon, estimator, truth, old_estimate = estimate, old_se = se)],
  result$results[, .(horizon, estimator, estimate, se, lower, upper)], by = c('horizon', 'estimator'))
fwrite(comparison, file.path(out, 'old-versus-expanded-library.csv'))
print(comparison, digits = 9)
print(result$ensemble_weights[kind %in% c('beta', 'adjoint')], digits = 7)
cat('Elapsed seconds:', result$seconds, '\n')
cat('Warnings:', paste(result$job_warnings, collapse = ' | '), '\n')
cat('Both SDR and TMLE completed. All weights retained. Original splits match exactly.\n')
