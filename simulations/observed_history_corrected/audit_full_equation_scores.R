# Separate read-only audit for full-equation cell-scored checks.
source(file.path(Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected"), "study.R"))
root <- Sys.getenv("SIM_OUTPUT", "results/observed-history-corrected/selection-revision")
scoring <- occupancy <- truncation <- support <- list()
all_points <- fread(file.path(root, "all-function-predictions.csv"))
sample_sizes <- as.integer(strsplit(Sys.getenv("SIM_SIZES", "4000,20000"), ",", fixed = TRUE)[[1L]])
for (mechanism in c("binary_longitudinal", "discrete_dose")) for (n in sample_sizes) {
  seed <- if (mechanism == "binary_longitudinal") 5103006L else 5203006L
  folder <- file.path(root, paste0(mechanism, "-n", n, "-seed", seed))
  set.seed(seed); g <- make_mechanism(mechanism, dose_max = 3L, visits = "missing_both")
  require_valid_mechanism(g)
  d <- draw_data(n, g)
  old <- if (n == 4000L) readRDS(file.path("results/observed-history-corrected/missing-both-study-v2/jobs",
    paste0(mechanism, "-n4000-r006.rds"))) else NULL
  result_file <- file.path(folder, "result.rds")
  result <- if (file.exists(result_file)) readRDS(result_file) else NULL
  reference <- Sys.getenv("SIM_REFERENCE")
  if (is.null(result) && nzchar(reference)) result <- readRDS(file.path(reference,
    paste0(mechanism, "-n", n, "-seed", seed), "result.rds"))
  prepared <- study_task(d, g, folds = if (!is.null(result)) result$folds else if (!is.null(old)) old$folds else 3L,
    learner_folds = if (!is.null(result)) result$learner_folds else if (!is.null(old)) old$learner_folds else NULL,
    learner_groups = 3L, balance_measurement = TRUE)
  task <- prepared$task; args <- study_arguments(g)
  for (fold in 1:3) {
    fit <- readRDS(file.path(folder, paste0("nuisance-fits-fold", fold, ".rds")))
    rows <- task$folds[[fold]]$training_set; valid <- task$folds[[fold]]$validation_set
    for (s in 1:2) {
      H <- task$vars$history("A", s); A <- task$vars$A[[s]]
      B <- as.matrix(task$natural[rows, c(H, A), drop = FALSE])
      V <- as.matrix(cbind(task$natural[rows, H, drop = FALSE], d[rows, args$health[[s]], drop = FALSE]))
      M <- d[[args$measurement[s]]][rows]
      measured_valid <- valid[d[[args$measurement[s]]][valid] == 1]
      newV <- as.matrix(cbind(task$natural[measured_valid, H, drop = FALSE],
        d[measured_valid, args$health[[s]], drop = FALSE]))
      raw <- predict(fit$beta_fits[[s]], newV)
      applied <- fit$beta_validation[d[[args$measurement[s]]][valid] == 1, s]
      stopifnot(max(abs(applied - pmin(100, pmax(-100, raw)))) < 1e-10)
      truncation[[length(truncation) + 1L]] <- data.table(mechanism, n, fold, horizon = s,
        measured_validation_n = length(raw), validation_clipped_n = sum(raw != applied),
        validation_max_change = max(abs(raw - applied)),
        training_values_at_limits = sum(abs(fit$beta_training[, s]) == 100, na.rm = TRUE))
      for (kind in c("beta", "adjoint")) {
        model <- fit[[paste0(kind, "_fits")]][[s]]
        stopifnot(kind != "adjoint" || !any(vapply(model$library, function(candidate)
          identical(candidate$method, "saturated_l1"), TRUE)))
        sieve <- model$candidates$sieve_md
        point_rows <- all_points$mechanism == mechanism & all_points$n == n &
          all_points$fold == fold & all_points$kind == kind &
          all_points$horizon == s & all_points$candidate == "sieve_md"
        target <- copy(all_points[point_rows])
        # A fitted dictionary includes only its own training combinations.
        target[, represented := cell %in% sieve$target_spec$levels]
        for (state in unique(target$current_R)) {
          ztarget <- target[current_R == state]
          squared <- ztarget$probability * (ztarget$estimate - ztarget$true)^2
          support[[length(support) + 1L]] <- data.table(mechanism, n, fold, kind,
            horizon = s, current_R = state, population_combinations = nrow(ztarget),
            represented_combinations = sum(ztarget$represented),
            missing_probability = sum(ztarget$probability[!ztarget$represented]) / sum(ztarget$probability),
            error_fraction_in_represented = if (sum(squared) > 0)
              sum(squared[ztarget$represented]) / sum(squared) else NA_real_)
        }
        Z <- if (kind == "beta") B else V
        old_v <- rep(0, ncol(model$cv_residuals))
        recomputed_u <- matrix(0, ncol(model$cv_residuals), ncol(model$cv_residuals))
        for (label in sort(unique(model$fold_id))) {
          index <- which(model$fold_id == label & if (kind == "beta") TRUE else M == 1)
          residual <- model$cv_residuals[index, , drop = FALSE]
          z <- Z[index, , drop = FALSE]
          metadata_score <- cmbridge::cell_moment_gram(residual, z)
          train_index <- which(model$fold_id != label & if (kind == "beta") TRUE else M == 1)
          spec <- cmbridge:::.ensemble_kernel(Z[train_index, , drop = FALSE], model$scoring_kernel$kernel,
            model$scoring_kernel$control, 1L + label)
          # Independent cellwise ordered-pair calculation, excluding self-products.
          stopifnot(model$scoring_kernel$kernel == "cell")
          groups <- base::split(seq_len(nrow(z)), cell_key(z))
          U <- Vscore <- matrix(0, ncol(residual), ncol(residual))
          for (group in groups) {
            r <- residual[group, , drop = FALSE]; k <- length(group)
            product <- tcrossprod(colSums(r))
            Vscore <- Vscore + product / (nrow(z) * k)
            if (k > 1L) U <- U + (product - crossprod(r)) / (nrow(z) * (k - 1L))
          }
          historical_v <- diag(Vscore)
          metadata <- attr(metadata_score, "cell_scoring")
          occupancy[[length(occupancy) + 1L]] <- data.table(mechanism, n, fold, kind, horizon = s,
            learner_fold = label, scored_n = length(index), cells = metadata$cells,
            singleton_cells = metadata$singleton_cells, paired_rows = metadata$paired_rows)
          weight <- length(index) / model$n_scored
          recomputed_u <- recomputed_u + weight * U
          # Historical V score for an audit only; never used to select a fit.
          old_v <- old_v + weight * historical_v
        }
        stopifnot(max(abs(recomputed_u - model$raw_gram)) < 1e-10,
          min(eigen(model$gram, symmetric = TRUE)$values) > -1e-10)
        scoring[[length(scoring) + 1L]] <- data.table(mechanism, n, fold, kind, horizon = s,
          candidate = names(model$weights), weight = as.numeric(model$weights),
          raw_u = diag(model$raw_gram), historical_v = old_v,
          removed_self_product_effect = old_v - diag(model$raw_gram),
          projected_loss = diag(model$gram), projection_norm = model$psd_projection$adjustment_norm)
      }
    }
  }
}
fwrite(rbindlist(scoring), file.path(root, "scorer-audit.csv"))
fwrite(rbindlist(occupancy), file.path(root, "scoring-cell-counts.csv"))
fwrite(rbindlist(truncation), file.path(root, "applied-bridge-truncation.csv"))
fwrite(rbindlist(support), file.path(root, "saturated-sieve-support.csv"))
cat("Verified raw U-statistics, PSD matrices, adjoint libraries, and applied bridge predictions.\n")
