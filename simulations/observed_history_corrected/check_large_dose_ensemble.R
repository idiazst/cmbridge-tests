# Reuse the original million-person sample and all its original assignments.
# True functions are evaluated only after the estimated-model fits finish.
source(file.path(Sys.getenv('STUDY_SOURCE'), 'study.R'))
out <- Sys.getenv('SIM_OUTPUT', 'results/observed-history-corrected/dose-ensemble-large')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path(out, 'large-fit.rds')
if (file.exists(checkpoint)) stop('Completed checkpoint already exists; preserve it.')
input <- 'results/observed-history-corrected/dose-cmbridge-debug/large-fit.rds'
reference <- readRDS(input)
stopifnot(reference$n == 1000000L, reference$seed == 5204001L, reference$fold == 2L,
          identical(names(cm_library()), c('sieve_md', 'landweber', 'pmmr', 'saturated_cv')))
g <- make_mechanism('discrete_dose', dose_max = 3L)
require_valid_mechanism(g)
set.seed(reference$seed); d <- draw_data(reference$n, g)
prepared <- study_task(d, g, folds = reference$folds, learner_folds = reference$learner_folds,
                       learner_groups = 3L, balance_measurement = TRUE)
task <- prepared$task; j <- reference$fold; rows <- task$folds[[j]]$training_set
stopifnot(identical(task$folds, reference$folds), identical(task$learner_folds, reference$learner_folds),
          length(rows) == reference$training_n)
args <- study_arguments(g); ctrl <- study_control(3L)
ctrl$.nested_split_function <- lmtp:::make_bridge_nested_folds(task$id[rows],
  d[rows, args$measurement, drop = FALSE], task$learner_folds[[j]])
original <- d; original$..bridge_row <- seq_len(nrow(d))
writeLines(c(capture.output(sessionInfo()),
  paste('original checkpoint SHA256', digest::digest(file = input, algo = 'sha256')),
  'Original sample: n=1000000, seed=5204001, outer training sample=2.',
  'All original estimator and learner assignments preserved.',
  'CM component check; ratio fits estimated; no sequential outcome regression run.'),
  file.path(out, 'session.txt'))
for (name in c('fit_bridge_ensemble', 'fit_adjoint_ensemble'))
  trace(name, where = asNamespace('cmbridge'), print = FALSE,
    tracer = quote(cat(format(Sys.time()), 'ensemble fit, training people:', length(M), '\n')))
trace('fit_cm', where = asNamespace('cmbridge'), print = FALSE,
  tracer = quote(cat(format(Sys.time()), 'candidate', method[1L], 'training rows:', length(response), '\n')))
started <- proc.time()[3L]; warnings <- character()
cat(format(Sys.time()), 'Starting original million-person nuisance fit.\n')
fitted <- withCallingHandlers(lmtp:::bridge_nuisance_fit(original, prepared$encoded, task, j,
  args$measurement, args$health, cm_library(FALSE, TRUE), cm_library(FALSE, FALSE),
  trt_library, TRUE, ctrl),
  warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
for (name in c('fit_bridge_ensemble', 'fit_adjoint_ensemble', 'fit_cm'))
  untrace(name, where = asNamespace('cmbridge'))
seconds <- proc.time()[3L] - started
cat(format(Sys.time()), 'Fits completed after', seconds, 'seconds. Evaluating exact population.\n')
p <- enumerate_data(g); pt <- study_task(p, g)$task
truth <- population_truth_functions(p, g)
ensemble <- predict_bridge_functions(fitted, pt, p, g)
H <- pt$vars$history('A', 2L); A <- pt$vars$A[[2L]]; M <- p$R3 == 1L
B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
V <- as.matrix(cbind(pt$natural[, H, drop = FALSE], p[, c('C3_covariate', 'Y3')]))
candidate <- lapply(seq_along(cm_library()), function(a) {
  beta <- rep(1, nrow(p)); beta[M] <- predict(fitted$beta_fits[[2L]]$candidates[[a]], V[M, , drop = FALSE])
  lambda <- predict(fitted$adjoint_fits[[2L]]$candidates[[a]], B)
  list(beta = beta, lambda = lambda)
})
names(candidate) <- names(cm_library())
stopifnot(all(vapply(candidate, function(x) all(is.finite(x$beta)) && all(is.finite(x$lambda)), TRUE)))
compact <- function(model) list(weights = model$weights, candidate_cv_loss = model$candidate_cv_loss,
  cv_loss = model$cv_loss, gram = model$gram, fold_gram = model$fold_gram,
  candidates = lapply(model$candidates, function(x) x[c('method', 'coefficients', 'target_spec',
    'instrument_spec', 'target_levels', 'instrument_levels', 'target_scale', 'instrument_scale',
    'target_nystrom', 'instrument_nystrom', 'penalty_cv', 'tuning', 'moment_loss', 'predict_fun')]))
weights <- rbindlist(lapply(c('beta', 'adjoint', 'ratio'), function(kind)
  fit_candidate_weights(fitted[[paste0(kind, '_fits')]], kind, j)), fill = TRUE)
counts <- data.table(cell = cell_key(task$natural[rows, c(task$vars$history('A', 2L),
  task$vars$A[[2L]]), drop = FALSE]), measured = d$R3[rows])[
    , .(training = .N, measured = sum(measured)), by = cell]
result <- list(n = reference$n, seed = reference$seed, fold = j, training_n = length(rows),
  beta = ensemble$beta, lambda = ensemble$lambda, candidate = candidate, population = p,
  beta_fit = compact(fitted$beta_fits[[2L]]), adjoint_fit = compact(fitted$adjoint_fits[[2L]]),
  folds = task$folds, learner_folds = task$learner_folds, counts = counts, ensemble_weights = weights,
  warnings = unique(warnings), seconds = seconds,
  cmbridge_version = as.character(packageVersion('cmbridge')), lmtp_version = as.character(packageVersion('lmtp')))
saveRDS(result, paste0(checkpoint, '.tmp')); stopifnot(file.rename(paste0(checkpoint, '.tmp'), checkpoint))
fwrite(weights, file.path(out, 'ensemble_weights.csv'))
fwrite(counts, file.path(out, 'training-counts.csv'))
cat('Bridge weights:\n'); print(result$beta_fit$weights, digits = 8)
cat('Adjoint weights:\n'); print(result$adjoint_fit$weights, digits = 8)
cat('Warnings:', paste(result$warnings, collapse = ' | '), '\n')
cat('Saved', checkpoint, '\n')
