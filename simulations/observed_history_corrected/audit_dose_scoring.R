# Diagnose and repair the large-sample score without refitting any learner.
library(data.table)
out <- 'results/observed-history-corrected/dose-cmbridge-debug'
inputs <- readRDS(file.path(out, 'scoring-inputs.rds'))
fit <- readRDS(file.path(out, 'large-fit.rds'))
candidate <- fit$beta_fit$candidates[[1L]]
index <- match(inputs$beta_key[inputs$observed], candidate$target_levels)
deviation <- numeric(length(index))
seen <- !is.na(index)
deviation[seen] <- candidate$coefficients[1L + index[seen]]
prediction <- pmax(candidate$tuning$lower, pmin(candidate$tuning$upper,
              candidate$coefficients[1L] + deviation))
residual <- rep(-1, length(inputs$observed))
residual[inputs$observed] <- prediction - 1

old_gram <- cmbridge:::.ensemble_gram
corrected_namespace <- new.env(parent = asNamespace('cmbridge'))
source('../cmbridge/R/ensemble.R', local = corrected_namespace)
new_gram <- corrected_namespace$.ensemble_gram
warning_record <- character()
old_score <- withCallingHandlers(old_gram(matrix(residual, ncol = 1L), inputs$B,
    list(kernel = 'cell')), warning = function(w) {
      warning_record <<- c(warning_record, conditionMessage(w),
                           paste(deparse(conditionCall(w)), collapse = ' '))
      invokeRestart('muffleWarning')
    })
new_score <- new_gram(matrix(residual, ncol = 1L), inputs$B, list(kernel = 'cell'))
groups <- match(cmbridge:::.cell_keys(inputs$B), unique(cmbridge:::.cell_keys(inputs$B)))
direct <- data.table(group = groups, residual = residual)[,
  .(rows = .N, mean_residual = mean(residual)), by = group][,
  sum(rows * mean_residual^2) / sum(rows)]
stopifnot(is.na(old_score[1L]), is.finite(new_score[1L]),
          abs(new_score[1L] - direct) < 1e-12,
          'n * tabulate(index)' %in% warning_record)

b_counts <- tabulate(groups)
v_counts <- as.integer(table(inputs$beta_key[inputs$observed]))
# All CM validation sets have at most one of the three root groups.
# In a two-group child fit, each side is also no larger than a root group.
max_validation_n <- max(tabulate(inputs$learner_labels))
maximum_cv_product_bound <- as.double(max_validation_n) * max(b_counts)
maximum_adjoint_product_bound <- as.double(sum(inputs$observed)) * max(v_counts)
stopifnot(maximum_cv_product_bound < .Machine$integer.max,
          maximum_adjoint_product_bound < .Machine$integer.max,
          all(is.finite(fit$tuning$loss)))
record <- data.table(
  old_final_bridge_training_score = old_score[1L],
  corrected_final_bridge_training_score = new_score[1L],
  direct_conditional_mean_score = direct,
  maximum_bridge_validation_count_product_bound = maximum_cv_product_bound,
  maximum_adjoint_count_product_bound = maximum_adjoint_product_bound,
  maximum_integer = .Machine$integer.max,
  overflowing_final_bridge_cells = sum(as.double(length(residual)) * b_counts > .Machine$integer.max),
  selected_penalty_losses_finite = all(is.finite(fit$tuning$loss)))
fwrite(record, file.path(out, 'integer-overflow-audit.csv'))
writeLines(warning_record, file.path(out, 'integer-overflow-warning.txt'))
print(record)
cat('Only the final bridge training-score arithmetic is affected. The code computes\n')
cat('that score after predictions, penalty selection, and weights have been determined.\n')
cat('The corrected score is obtained from the saved fitted coefficients; no refit is required.\n')

# Summaries of the remaining error in the large sample.
bp <- fread(file.path(out, 'bridge-predictions.csv'))[source == 'large check']
ap <- fread(file.path(out, 'adjoint-equation-predictions.csv'))[source == 'large check']
for (kind in c('bridge', 'adjoint-equation')) {
  points <- if (kind == 'bridge') copy(bp) else copy(ap)
  points[, probability := probability / sum(probability)]
  points[, `:=`(error = estimate - true,
                squared_error_contribution = probability * (estimate - true)^2)]
  points[, training_measured := as.integer(table(inputs$beta_key[inputs$observed])[cell])]
  stopifnot(all(!is.na(points$training_measured)))
  fwrite(points[order(-abs(error))], file.path(out, paste0(kind, '-large-errors.csv')))
}
