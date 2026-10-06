# Repair a metadata assignment in the completed run. No fit or prediction changes.
library(data.table)
out <- 'results/observed-history-corrected/dose-ensemble-r002'
path <- file.path(out, 'discrete_dose-n4000-r002.rds')
new <- readRDS(path)
original <- file.path(out, 'discrete_dose-n4000-r002-executed.rds')
if (!file.exists(original)) stopifnot(file.copy(path, original))
weights <- copy(new$ensemble_weights)
rows <- which(weights$kind == 'regression')
stopifnot(length(rows) %% 2L == 0L)
blocks <- split(rows, rep(c('sdr', 'tmle'), each = length(rows) / 2L))
for (name in names(blocks)) {
  ids <- blocks[[name]]
  stopifnot(identical(unique(weights$fold[ids]), 1:3))
  set(weights, i = ids, j = 'estimator', value = name)
}
stopifnot(all(weights[, abs(sum(weight) - 1) < 1e-8,
  by = .(kind, fold, horizon_or_depth, estimator)]$V1))
# Validate every repaired label against the independently saved saturated
# candidate tuning records, whose estimator labels were already correct.
check <- merge(weights[kind == 'regression' & method == 'SL.saturated_l1',
    .(fold, horizon_or_depth, estimator, saved_weight = weight)],
  unique(new$tuning[kind == 'regression' & selected == TRUE,
    .(fold, horizon_or_depth, estimator, tuning_weight = weight)]),
  by = c('fold', 'horizon_or_depth', 'estimator'))
stopifnot(nrow(check) == length(rows) / 3L,
          max(abs(check$saved_weight - check$tuning_weight)) < 1e-12)
new$ensemble_weights <- weights
stopifnot(identical(new$results, readRDS(original)$results),
          identical(new$fitted_functions, readRDS(original)$fitted_functions))
saveRDS(new, path)
fwrite(weights, file.path(out, 'ensemble_weights.csv'))
cat('Repaired only SDR/TMLE labels in regression weight metadata.\n')
cat('All labels independently checked against the saved penalty-tuning records.\n')
cat('Estimates, fitted functions, candidate weights, and validation losses unchanged.\n')
