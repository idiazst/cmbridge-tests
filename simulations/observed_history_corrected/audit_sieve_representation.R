# Algebraic checks of the saved sieve classes and empirical moment operators.
# No refitting, new datasets, zero-penalty fits, or known-function simulations.
source(file.path(Sys.getenv("STUDY_SOURCE"), "study.R"))
root <- Sys.getenv("SIM_OUTPUT")
points <- fread(file.path(root, "all-function-predictions.csv"))
classes <- blocks <- list()
for (mechanism in c("binary_longitudinal", "discrete_dose")) for (sample_n in c(4000L, 20000L)) {
  seed <- if (mechanism == "binary_longitudinal") 5103006L else 5203006L
  folder <- file.path(root, paste0(mechanism, "-n", sample_n, "-seed", seed))
  result <- readRDS(file.path(folder, "result.rds"))
  g <- make_mechanism(mechanism, dose_max = 3L, visits = "missing_both")
  task <- study_task(result$data, g, folds = result$folds,
    learner_folds = result$learner_folds, learner_groups = 3L)$task
  arguments <- study_arguments(g)
  for (fold in 1:3) {
    fits <- readRDS(file.path(folder, paste0("nuisance-fits-fold", fold, ".rds")))
    training <- task$folds[[fold]]$training_set
    for (s in 1:2) {
      fit <- fits$beta_fits[[s]]$candidates$sieve_md
      selected <- points$mechanism == mechanism & points$n == sample_n &
        points$fold == fold & points$kind == "beta" & points$horizon == s &
        points$candidate == "sieve_md"
      z <- copy(points[selected])
      truth_seen <- z$true[match(fit$target_spec$levels, z$cell)]
      stopifnot(!anyNA(truth_seen), all(truth_seen > 1))
      # Reuse the actual fitted prediction code with algebraically chosen
      # coefficients. This checks representability, not an oracle estimator.
      f <- fit$predict_fun
      environment(f) <- list2env(list(xspec = fit$target_spec,
        link = "inverse_logit", link_coef = -log(truth_seen - 1),
        link_intercept = -log(mean(truth_seen) - 1)),
        parent = parent.env(environment(f)))
      x <- do.call(rbind, lapply(strsplit(z$cell, "|", fixed = TRUE), as.numeric))
      represented <- z$cell %in% fit$target_spec$levels
      reproduced <- f(x)
      stopifnot(max(abs(reproduced[represented] - z$true[represented])) < 1e-12)
      for (state in unique(z$current_R)) {
        ii <- z$current_R == state
        classes[[length(classes) + 1L]] <- data.table(mechanism, n = sample_n,
          fold, horizon = s, current_R = state,
          population_combinations = sum(ii), represented_combinations = sum(ii & represented),
          maximum_representation_error_seen = max(abs(reproduced[ii & represented] - z$true[ii & represented])),
          maximum_representation_error_full = max(abs(reproduced[ii] - z$true[ii])),
          omitted_probability = sum(z$probability[ii & !represented]) / sum(z$probability[ii]),
          lambda = fit$tuning$lambda, weight_ridge = fit$tuning$weight_ridge)
      }
      H <- task$vars$history("A", s); A <- task$vars$A[[s]]
      B <- as.matrix(task$natural[training, c(H, A), drop = FALSE])
      V <- as.matrix(cbind(task$natural[training, H, drop = FALSE],
        result$data[training, arguments$health[[s]], drop = FALSE]))
      active <- result$data[[arguments$measurement[s]]][training] == 1
      bkey <- cell_key(B); vkey <- rep(NA_character_, length(training))
      vkey[active] <- cell_key(V[active, , drop = FALSE])
      hkey <- cell_key(B[, seq_along(H), drop = FALSE])
      population_h <- cell_key(x[, seq_along(H), drop = FALSE])
      for (h in unique(population_h)) {
        ii <- which(hkey == h); jj <- which(population_h == h)
        target_levels <- unique(vkey[ii[active[ii]]])
        target_levels <- target_levels[!is.na(target_levels)]
        instrument_levels <- unique(bkey[ii])
        measured <- ii[active[ii]]
        rank <- 0L; condition <- Inf; smallest <- NA_real_
        if (length(target_levels) && length(instrument_levels)) {
          J <- as.matrix(Matrix::sparseMatrix(i = match(bkey[measured], instrument_levels),
            j = match(vkey[measured], target_levels), x = 1,
            dims = c(length(instrument_levels), length(target_levels))))
          counts <- tabulate(match(bkey[ii], instrument_levels), nbins = length(instrument_levels))
          singular <- svd(J / sqrt(counts), nu = 0, nv = 0)$d
          rank <- sum(singular > max(singular) * 1e-10)
          if (rank == length(target_levels)) {
            smallest <- min(singular); condition <- max(singular) / smallest
          }
        }
        blocks[[length(blocks) + 1L]] <- data.table(mechanism, n = sample_n,
          fold, horizon = s, current_R = z$current_R[jj[1L]], history = h,
          training_n = length(ii), measured_training_n = length(measured),
          represented_targets = length(target_levels), instrument_combinations = length(instrument_levels),
          operator_rank = rank, rank_deficient = rank < length(target_levels),
          all_population_targets_represented = length(target_levels) == length(jj),
          condition_number = condition, smallest_singular_value = smallest,
          population_probability = sum(z$probability[jj]),
          fitted_squared_error = sum(z$probability[jj] * (z$estimate[jj] - z$true[jj])^2))
      }
    }
  }
}
fwrite(rbindlist(classes), file.path(root, "sieve-representation-audit.csv"))
fwrite(rbindlist(blocks), file.path(root, "sieve-empirical-operator-audit.csv"))
cat("The actual inverse-expit prediction code reproduces the truth on every represented bridge combination to 1e-12.\n")
cat("Saved training dictionaries omit population combinations; inspect the full-support errors separately.\n")
