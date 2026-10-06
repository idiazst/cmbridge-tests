source(file.path(Sys.getenv('STUDY_SOURCE', 'simulations/observed_history_corrected'), 'study.R'))
out <- 'results/observed-history-corrected/dose-ensemble-r002'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
old <- readRDS('results/observed-history-corrected/dose-seed-check-r002/jobs/discrete_dose-n4000-r002.rds')
set.seed(old$design$seed)
g <- make_mechanism('discrete_dose', dose_max = 3L)
d <- draw_data(old$design$n, g)
prepared <- study_task(d, g, folds = old$folds, learner_folds = old$learner_folds,
                       learner_groups = 3L, balance_measurement = TRUE)
task <- prepared$task; rows <- old$folds[[2L]]$training_set
H <- task$vars$history('A', 2L); A <- task$vars$A[[2L]]
B <- as.matrix(task$natural[rows, c(H, A), drop = FALSE])
V <- as.matrix(cbind(task$natural[rows, H, drop = FALSE], d[rows, c('C3_covariate', 'Y3')]))
M <- d$R3[rows]
control <- list(target_basis = 'cell', instrument_basis = 'cell', lower = 1, upper = 6)
records <- list()
for (iteration in seq_len(3L)) for (method in c('sieve_md', 'landweber')) {
  elapsed <- system.time(model <- cmbridge::fit_bridge(B, V, M, method, control))[['elapsed']]
  records[[length(records) + 1L]] <- data.table(method = method, iteration = iteration,
    seconds = elapsed, training_n = length(rows), target_cells = length(model$target_spec$levels),
    conditioning_cells = length(model$instrument_spec$levels))
}
records <- rbindlist(records)
summary <- records[, .(median_seconds = median(seconds)), by = method][order(median_seconds)]
fwrite(records, file.path(out, 'joint-basis-timing.csv'))
fwrite(summary, file.path(out, 'joint-basis-timing-summary.csv'))
print(summary)
cat('Faster method:', summary$method[1L], '\n')
