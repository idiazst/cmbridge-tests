# Shared implementation for the approved, estimated-model study. Exact
# population calculations below are diagnostics only, never fitting inputs.
library(data.table)
library(lmtp)
future::plan(future::sequential)
data.table::setDTthreads(1L)
stopifnot(packageVersion("lmtp") == "1.6.0.9021",
          packageVersion("cmbridge") == "0.3.0.9021")
if (!requireNamespace("earth", quietly = TRUE))
  stop("The regression library requires the earth package for MARS.")
study_source <- Sys.getenv("STUDY_SOURCE", "simulations/observed_history_corrected")
source(file.path(study_source, "dgp.R"))

study_control <- function(learner_groups = 2L) {
  control <- lmtp_control(.learners_trt_folds = learner_groups,
  .learners_outcome_folds = learner_groups, .discrete = FALSE, .trim = 1,
  .return_full_fits = TRUE, .isotonic_constraint = FALSE, .cm_kernel = "rbf")
  control$.bridge_prediction_bounds <- c(-100, 100)
  control
}
study_shift <- function(g) function(data, trt) {
  out <- data[, trt, drop = FALSE]
  time <- as.integer(sub("A([0-9]+)_.*", "\\1", trt[1]))
  measured <- data[[paste0("R", time)]] == 1
  out[[trt[1]]][measured] <- if (g$dose) pmin(out[[trt[1]]][measured] + 1, max(g$A$policy_component)) else 1
  out
}
study_arguments <- function(g) list(
  trt = lapply(seq_len(g$tau), function(t) paste0("A", t, c("_policy", "_other"))),
  outcome = paste0("Y", 2:(g$tau + 1L)),
  measurement = paste0("R", 2:(g$tau + 1L)),
  health = lapply(2:(g$tau + 1L), function(t) c(paste0("C", t, "_covariate"), paste0("Y", t))),
  baseline = "C1_baseline",
  time_vary = lapply(seq_len(g$tau), function(t) c(paste0("C", t, "_covariate"), paste0("Y", t), paste0("R", t))),
  shift = study_shift(g), mtp = g$dose, bounds = c(-100, 100), control = study_control())
study_task <- function(d, g, folds = 3, learner_folds = NULL, learner_groups = 2L, balance_measurement = FALSE) {
  args <- study_arguments(g); args$control <- study_control(learner_groups); e <- as.data.frame(d)
  e$..bridge_row <- seq_len(nrow(e)); L <- list()
  for (t in seq_len(g$tau)) {
    L[[t]] <- character()
    for (name in args$time_vary[[t]]) {
      key <- paste0("..bridge_", t, "_", name); miss <- paste0(key, "_missing")
      e[[key]] <- ifelse(is.na(d[[name]]), 0, d[[name]])
      e[[miss]] <- as.numeric(is.na(d[[name]]))
      L[[t]] <- c(L[[t]], key, miss)
    }
  }
  task <- lmtp:::LmtpTask$new(e, lmtp:::make_shifted(e, args$trt, NULL, args$shift, NULL),
    args$trt, args$outcome, L, args$baseline, NULL, NULL, Inf, NULL,
    "continuous", folds, NULL, c(-100, 100), learner_folds)
  lmtp:::prepare_curve_splits(task, args$control)
  if (balance_measurement && is.null(learner_folds)) for (j in seq_along(task$folds)) {
    rows <- task$folds[[j]]$training_set
    profile <- do.call(paste, c(d[rows, args$measurement, drop = FALSE], sep = "|"))
    ids <- as.character(task$id[rows]); labels <- integer(length(rows)); offset <- 0L
    for (stratum in sort(unique(profile))) {
      group <- which(profile == stratum); group <- group[order(ids[group])]
      labels[group] <- (seq_along(group) - 1L + offset) %% learner_groups + 1L
      offset <- offset + length(group)
    }
    task$learner_folds[[j]] <- labels
  }
  list(task = task, encoded = e)
}
cm_library <- function(constant = FALSE, bridge = TRUE) {
  if (constant) stop("Misspecified conditional-moment libraries are deferred in this revision.")
  library <- list(sieve_md = list(method = "sieve_md",
    control = list(target_basis = if (bridge) "cell_linear" else "cell",
      instrument_basis = if (bridge) "poly" else "cell", instrument_degree = 1L,
      penalty_on = if (bridge) "link" else "function",
      link = if (bridge) "inverse_logit" else "identity",
      penalty_scales = 10^seq(-10, 0, by = .25))),
    landweber = list(method = "landweber",
      control = if (bridge) list(link = "inverse_logit", target_basis = "poly", target_degree = 1L,
        instrument_basis = "poly", instrument_degree = 1L, precondition = TRUE) else list(link = "identity")),
    pmmr = list(method = "pmmr",
      control = list(link = if (bridge) "inverse_logit" else "identity")))
  library
}
sl_library <- c("SL.saturated_l1_cv", "SL.earth", "SL.mean")
trt_library <- c("SL.saturated_l1_cv", "SL.mean")

predict_ratios <- function(fits, task, g) {
  n <- nrow(task$natural); answer <- matrix(0, n, g$tau)
  for (t in seq_len(g$tau)) {
    A <- task$vars$A[[t]]
    pred <- predict(fits[[t]], task$natural)
    followed <- lmtp:::followed_rule(task$natural, task$shifted, A, g$dose)
    if (!g$dose) pred[followed] <- pmax(pred[followed], .5)
    answer[, t] <- pred * followed / (1 - pmin(pred, .999))
  }
  answer
}
fit_ratios_without_history <- function(task, fold, g, split_cache, split_function) {
  rows <- task$folds[[fold]]$training_set; valid <- task$folds[[fold]]$validation_set
  labels <- task$learner_folds[[fold]]
  fit_on <- function(train, validation) {
    fits <- vector("list", g$tau); out <- matrix(0, length(validation), g$tau)
    for (t in seq_len(g$tau)) {
      A <- task$vars$A[[t]]
      stacked <- lmtp:::stack_data(task$natural[rows[train], ], task$shifted[rows[train], ],
                                  task$vars$A, task$vars$C, t)
      # Only current treatment is retained in this deliberately wrong model.
      vars <- c("..i..lmtp_id", A, "..i..lmtp_stack_indicator")
      ids <- stacked$..i..lmtp_id
      if (length(unique(labels[train])) > 1L) inner <- labels[train] else {
        key <- digest::digest(train)
        if (!exists(key, split_cache, inherits = FALSE)) {
          plan <- lmtp:::make_folds(task$natural[rows[train], ], V = 2, cluster_ids = task$id[rows[train]])
          inner <- integer(length(train))
          for (k in seq_along(plan)) inner[plan[[k]]$validation_set] <- k
          assign(key, inner, split_cache)
        }
        inner <- get(key, split_cache, inherits = FALSE)
      }
      inner <- inner[match(ids, task$id[rows[train]])]
      fits[[t]] <- lmtp:::run_ensemble(stacked[, vars], "..i..lmtp_stack_indicator",
        trt_library, "binomial", "..i..lmtp_id", list(fold_id = inner), FALSE, FALSE,
        penalty_folds = split_function)
      pred <- predict(fits[[t]], task$natural[validation, ])
      followed <- lmtp:::followed_rule(task$natural[validation, ], task$shifted[validation, ], A, g$dose)
      if (!g$dose) pred[followed] <- pmax(pred[followed], .5)
      out[, t] <- pred * followed / (1 - pmin(pred, .999))
    }
    list(ratios = out, fits = fits)
  }
  final <- fit_on(seq_along(rows), valid)
  training <- matrix(NA_real_, length(rows), g$tau)
  for (label in unique(labels)) {
    vi <- which(labels == label)
    training[vi, ] <- fit_on(which(labels != label), rows[vi])$ratios
  }
  list(ratios_training = training, ratios_validation = final$ratios, ratio_fits = final$fits)
}

state_index <- function(d, t) ifelse(d[[paste0("R", t)]] == 0, 1L,
  2L + ifelse(is.na(d[[paste0("C", t, "_covariate")]]), 0, d[[paste0("C", t, "_covariate")]]) +
    2L * ifelse(is.na(d[[paste0("Y", t)]]), 0, d[[paste0("Y", t)]]))
action_index <- function(d, t, g) 1L + d[[paste0("A", t, "_policy")]] +
  (max(g$A$policy_component) + 1L) * d[[paste0("A", t, "_other")]]
enumerate_data <- function(g) {
  d <- data.frame(C1_baseline = 0:1, C1_covariate = 0, Y1 = 0, R1 = 1, probability = .5)
  if (g$visits == "missing_both") d$Y1 <- d$C1_baseline
  for (t in seq_len(g$tau)) {
    pieces <- list()
    for (a in seq_len(nrow(g$A))) for (c in seq_len(nrow(g$C))) for (m in 0:1) {
      x <- d; h <- state_index(d, t)
      probability <- g$g[cbind(d$C1_baseline + 1L, t, h, a)] *
        g$health[cbind(d$C1_baseline + 1L, t, h, a, c)] *
        dbinom(m, 1, g$measurement[cbind(d$C1_baseline + 1L, t, h, c)])
      x$probability <- d$probability * probability
      x[[paste0("A", t, "_policy")]] <- g$A$policy_component[a]
      x[[paste0("A", t, "_other")]] <- g$A$other_component[a]
      x[[paste0("R", t + 1L)]] <- m
      x[[paste0("C", t + 1L, "_covariate")]] <- if (m) g$C$covariate[c] else NA_real_
      x[[paste0("Y", t + 1L)]] <- if (m) g$C$outcome[c] else NA_real_
      pieces[[length(pieces) + 1L]] <- x[x$probability > 0, ]
    }
    d <- as.data.table(do.call(rbind, pieces))
    vars <- setdiff(names(d), "probability")
    d <- as.data.frame(d[, .(probability = sum(probability)), by = vars])
  }
  stopifnot(abs(sum(d$probability) - 1) < 1e-12)
  d
}
population_truth_functions <- function(p, g) {
  n <- nrow(p); beta <- ratios <- lambda <- matrix(1, n, g$tau); previous <- rep(1, n)
  for (t in seq_len(g$tau)) {
    h <- state_index(p, t); a <- action_index(p, t, g)
    ratios[, t] <- g$ratio[cbind(p$C1_baseline + 1L, t, h, a)]
    measured <- p[[paste0("R", t + 1L)]] == 1
    c <- 1L + p[[paste0("C", t + 1L, "_covariate")]][measured] +
      2L * p[[paste0("Y", t + 1L)]][measured]
    beta[measured, t] <- g$beta[cbind(p$C1_baseline[measured] + 1L, t, h[measured], c)]
    for (b in 1:2) for (state in unique(h[p$C1_baseline == b - 1L])) {
      joint <- g$health[b, t, state, , ] * g$g[b, t, state, ]
      # Some next-health values have zero probability at a missed visit.
      # Conditional equations apply only to the reachable support.
      active <- which(colSums(joint) > 0)
      posterior <- t(joint[, active, drop = FALSE]) / colSums(joint)[active]
      target <- g$C$outcome[active] * as.vector(posterior %*% g$ratio[b, t, state, ])
      sv <- svd(posterior); keep <- sv$d > 1e-10 * max(sv$d)
      ell <- as.vector(sv$v[, keep, drop = FALSE] %*%
        (as.vector(crossprod(sv$u[, keep, drop = FALSE], target)) / sv$d[keep]))
      ii <- which(p$C1_baseline == b - 1L & h == state)
      lambda[ii, t] <- previous[ii] * ell[a[ii]]
    }
    previous <- previous * ratios[, t]
  }
  list(beta = beta, lambda = lambda, ratios = ratios)
}
cell_key <- function(x) cmbridge:::.cell_keys(as.matrix(x))
conditional_mean <- function(y, w, key) {
  value <- data.table(cell = key, y = y, w = w)[, .(value = sum(w * y) / sum(w)), by = cell]
  setNames(value$value, value$cell)
}
raw_scores <- function(d, g, nu) {
  score <- sequential <- matrix(0, nrow(d), g$tau)
  for (s in seq_len(g$tau)) {
    measured <- d[[paste0("R", s + 1L)]] == 1
    terminal <- numeric(nrow(d))
    terminal[measured] <- d[[paste0("Y", s + 1L)]][measured] * nu$beta[measured, s]
    omega <- rep(1, nrow(d)); value <- nu$shifted[[s]][, 1]
    for (t in seq_len(s)) {
      omega <- omega * nu$ratios[, t]
      after <- if (t == s) terminal else nu$shifted[[s]][, t + 1L]
      value <- value + omega * (after - nu$natural[[s]][, t])
    }
    sequential[, s] <- value
    residual <- rep(-1, nrow(d)); residual[measured] <- nu$beta[measured, s] - 1
    score[, s] <- value - nu$lambda[, s] * residual
  }
  list(score = score, sequential = sequential)
}
predict_bridge_functions <- function(fitted, pt, p, g) {
  beta <- matrix(1, nrow(p), g$tau); lambda <- matrix(0, nrow(p), g$tau)
  for (s in seq_len(g$tau)) {
    if (is.null(fitted$beta_fits[[s]])) next
    H <- pt$vars$history("A", s); A <- pt$vars$A[[s]]
    B <- as.matrix(pt$natural[, c(H, A), drop = FALSE])
    V <- as.matrix(cbind(pt$natural[, H, drop = FALSE],
      p[, c(paste0("C", s + 1L, "_covariate"), paste0("Y", s + 1L)), drop = FALSE]))
    observed <- p[[paste0("R", s + 1L)]] == 1
    beta[observed, s] <- lmtp:::bridge_predict(fitted$beta_fits[[s]], V[observed, , drop = FALSE], study_control())
    lambda[, s] <- predict(fitted$adjoint_fits[[s]], B)
  }
  list(beta = beta, lambda = lambda)
}
predict_regressions <- function(fits, task, g, method, targeting = NULL) {
  natural <- lmtp:::pivot(task$natural, task$vars)
  shifted <- natural; other <- lmtp:::pivot(task$shifted, task$vars)
  A <- grep("^..i..A", names(shifted), value = TRUE)
  A <- grep("_lag", A, value = TRUE, invert = TRUE)
  shifted[, A] <- other[, A]
  qn <- lapply(seq_len(g$tau), function(s) matrix(0, task$n, s))
  qd <- lapply(seq_len(g$tau), function(s) matrix(0, task$n, s + 1L))
  for (depth in seq_len(g$tau)) {
    tol <- if (method == "tmle") study_control()$.bound else NULL
    qn <- lmtp:::update_m(qn, fits[[depth]], natural, depth, g$tau, FALSE, tol)
    qd <- lmtp:::update_m(qd, fits[[depth]], shifted, depth, g$tau, FALSE, tol)
  }
  if (method == "tmle") for (s in seq_len(g$tau)) for (r in rev(seq_len(s))) {
    epsilon <- targeting[[s]]$epsilon[targeting[[s]]$time == r]
    bound <- function(x) pmin(1 - study_control()$.bound, pmax(study_control()$.bound, x))
    qn[[s]][, r] <- plogis(qlogis(bound(qn[[s]][, r])) + epsilon)
    qd[[s]][, r] <- plogis(qlogis(bound(qd[[s]][, r])) + epsilon)
  }
  list(natural = lapply(qn, task$rescale),
       shifted = lapply(seq_len(g$tau), function(s) task$rescale(qd[[s]][, seq_len(s), drop = FALSE])))
}
population_diagnostics <- function(p, pt, g, nu, truth_nu, truth) {
  w <- p$probability; scores <- raw_scores(p, g, nu); ans <- list()
  for (s in seq_len(g$tau)) {
    measured <- p[[paste0("R", s + 1L)]] == 1
    terminal <- numeric(nrow(p)); terminal[measured] <- p[[paste0("Y", s + 1L)]][measured] * nu$beta[measured, s]
    omega <- apply(truth_nu$ratios[, seq_len(s), drop = FALSE], 1, prod)
    beta_target <- sum(w * omega * terminal)
    seq_error <- sum(w * scores$sequential[, s]) - beta_target
    bridge_error <- sum(w * scores$score[, s]) - truth[s] - seq_error
    product <- -sum(w[measured] * (nu$lambda[measured, s] - truth_nu$lambda[measured, s]) *
      (nu$beta[measured, s] - truth_nu$beta[measured, s]))
    stopifnot(abs(bridge_error - product) < 1e-8)
    current <- terminal; ratio_error <- m_error <- numeric(s)
    for (t in s:1) {
      H <- pt$vars$history("A", t); A <- pt$vars$A[[t]]
      B <- pt$natural[, c(H, A), drop = FALSE]
      means <- conditional_mean(current, w, cell_key(B))
      exact <- unname(means[cell_key(B)])
      m_error[t] <- sum(w * (nu$natural[[s]][, t] - exact)^2)
      ratio_error[t] <- sum(w * (nu$ratios[, t] - truth_nu$ratios[, t])^2)
      B[, A] <- pt$shifted[, A, drop = FALSE]
      current <- unname(means[cell_key(B)])
    }
    H <- pt$vars$history("A", s); A <- pt$vars$A[[s]]
    B <- pt$natural[, c(H, A), drop = FALSE]
    V <- cbind(pt$natural[, H, drop = FALSE], p[, c(paste0("C", s + 1L, "_covariate"), paste0("Y", s + 1L))])
    residual <- rep(-1, nrow(p)); residual[measured] <- nu$beta[measured, s] - 1
    bridge_moments <- conditional_mean(residual, w, cell_key(B))
    omega_y <- omega * ifelse(measured, p[[paste0("Y", s + 1L)]], 0)
    adjoint_moments <- conditional_mean(nu$lambda[measured, s] - omega_y[measured],
      w[measured], cell_key(V[measured, , drop = FALSE]))
    ans[[s]] <- data.frame(horizon = s, population_bias = sum(w * scores$score[, s]) - truth[s],
      sequential_remainder = seq_error, bridge_remainder = bridge_error,
      bridge_identity_error = bridge_error - product,
      sequential_product_norm = sum(sqrt(ratio_error * m_error)),
      beta_mse_measured = sum(w[measured] * (nu$beta[measured, s] - truth_nu$beta[measured, s])^2) / sum(w[measured]),
      lambda_mse_to_selected_solution = sum(w * (nu$lambda[, s] - truth_nu$lambda[, s])^2),
      ratio_mse_sum = sum(ratio_error), regression_mse_sum = sum(m_error),
      bridge_equation_max_error = max(abs(bridge_moments)),
      adjoint_equation_max_error = max(abs(adjoint_moments)),
      beta_lower_fraction = sum(w[measured] * (nu$beta[measured, s] <= 1 + 1e-6)) / sum(w[measured]),
      beta_upper_fraction = sum(w[measured] * (nu$beta[measured, s] >= 6 - 1e-6)) / sum(w[measured]))
  }
  rbindlist(ans)
}
fit_tuning <- function(fitted, kind, fold) {
  ans <- list()
  for (s in seq_along(fitted)) if (!is.null(fitted[[s]])) {
    object <- fitted[[s]]
    if (inherits(object, "lmtp_ensemble")) {
      if (any(object$errorsInCVLibrary) || any(object$errorsInLibrary))
        stop("A SuperLearner candidate failed; this comparison cannot count as a successful fit.")
      candidates <- object$fitLibrary; weights <- object$coef
    } else {candidates <- object$candidates; weights <- object$weights}
    for (a in seq_along(candidates)) {
      cv <- candidates[[a]]$penalty_cv
      if (is.null(cv)) next
      cv <- as.data.table(cv)
      cv[, `:=`(kind = kind, fold = fold, horizon_or_depth = s, candidate = a, weight = weights[a],
        penalty_status = candidates[[a]]$penalty_status,
        boundary_failures = nrow(candidates[[a]]$penalty_boundary))]
      ans[[length(ans) + 1L]] <- cv
    }
  }
  rbindlist(ans, fill = TRUE)
}
fit_candidate_weights <- function(fitted, kind, fold) {
  ans <- list()
  for (s in seq_along(fitted)) if (!is.null(fitted[[s]])) {
    object <- fitted[[s]]
    if (inherits(object, "lmtp_ensemble")) {
      if (any(object$errorsInCVLibrary) || any(object$errorsInLibrary))
        stop("A SuperLearner candidate failed; cannot report successful ensemble validation.")
      weights <- object$coef; loss <- object$cvRisk
      candidate <- names(weights)
      method <- vapply(object$fitLibrary, function(x) class(x)[1L], "")
    } else {
      weights <- object$weights; loss <- object$candidate_cv_loss
      candidate <- names(weights)
      method <- vapply(object$library, function(x) x$method, "")
    }
    ans[[length(ans) + 1L]] <- data.table(kind = kind, fold = fold,
      horizon_or_depth = s, candidate = candidate, method = unname(method),
      weight = as.numeric(weights), cv_loss = as.numeric(loss))
  }
  rbindlist(ans, fill = TRUE)
}
joint_critical <- function(scores, draws = 2000L) {
  covariance <- cov(scores)
  if (ncol(scores) == 1L) covariance <- matrix(covariance, 1, 1)
  variance <- diag(covariance)
  if (any(!is.finite(variance)) || any(variance <= 0)) stop("Nonpositive gradient variance.")
  correlation <- covariance / sqrt(outer(variance, variance))
  eig <- eigen(correlation, symmetric = TRUE)
  z <- matrix(rnorm(draws * ncol(scores)), draws, ncol(scores)) %*%
    (diag(sqrt(pmax(0, eig$values)), ncol(scores)) %*% t(eig$vectors))
  z <- z * sqrt((nrow(scores) - 1) / nrow(scores))
  # This samples the exact conditional law of Gaussian person multipliers
  # contracted with the empirical gradient, avoiding an n x 2000 matrix.
  as.numeric(quantile(apply(abs(z), 1, max), .95, names = FALSE))
}

run_dataset <- function(mechanism, n, replicate, seed, all_correct = FALSE, dose_max = 2L, learner_groups = 2L, balance_measurement = FALSE,
                        folds = 3L, learner_folds = NULL, visits = "scheduled", store_scores = FALSE, nuisance_callback = NULL, nuisance_cache = NULL) {
  set.seed(seed); g <- make_mechanism(mechanism, dose_max = dose_max, visits = visits); require_valid_mechanism(g)
  truth <- exact_truth(g); stopifnot(max(abs(truth - backward_truth(g))) < 1e-12)
  d <- draw_data(n, g); prepared <- study_task(d, g, folds = folds, learner_folds = learner_folds,
    learner_groups = learner_groups, balance_measurement = balance_measurement); task <- prepared$task
  p <- enumerate_data(g); pt <- study_task(p, g)$task
  truth_nu <- population_truth_functions(p, g)
  args <- study_arguments(g); args$control <- study_control(learner_groups); original <- d; original$..bridge_row <- seq_len(n)
  fitted <- population <- wrong_ratio <- vector("list", length(task$folds))
  tuning <- occupancy <- ensemble_weights <- list()
  for (j in seq_along(task$folds)) {
    shared_ratios <- new.env(parent = emptyenv()); shared_splits <- new.env(parent = emptyenv())
    train <- task$folds[[j]]$training_set
    if (is.function(nuisance_callback)) cat(format(Sys.time()), "Fitting nuisance functions for outer split", j, "\n")
    nested <- lmtp:::make_bridge_nested_folds(task$id[train], d[train, args$measurement, drop = FALSE],
                                              task$learner_folds[[j]])
    fold_control <- args$control; fold_control$.nested_split_function <- nested
    fitted[[j]] <- population[[j]] <- vector("list", 2)
    for (constant in if (all_correct) FALSE else c(FALSE, TRUE)) {
      k <- as.integer(constant) + 1L
      fitted[[j]][[k]] <- if (is.function(nuisance_cache)) nuisance_cache(j, constant) else NULL
      if (is.null(fitted[[j]][[k]])) fitted[[j]][[k]] <- lmtp:::bridge_nuisance_fit(original, prepared$encoded, task, j,
        args$measurement, args$health, cm_library(constant, TRUE), cm_library(constant, FALSE),
        trt_library, g$dose, fold_control, shared_ratios, shared_splits)
      fitted[[j]][[k]]$response_nuisance_cache <- new.env(parent = emptyenv())
      population[[j]][[k]] <- predict_bridge_functions(fitted[[j]][[k]], pt, p, g)
      population[[j]][[k]]$ratios <- predict_ratios(fitted[[j]][[k]]$ratio_fits, pt, g)
      if (is.function(nuisance_callback)) nuisance_callback(fitted[[j]][[k]], task, p, pt, g, j)
      for (kind in c("beta", "adjoint", "ratio")) {
        objects <- fitted[[j]][[k]][[paste0(kind, "_fits")]]
        tt <- fit_tuning(objects, kind, j)
        if (nrow(tt)) {tt[, constant := constant]; tuning[[length(tuning) + 1L]] <- tt}
        ww <- fit_candidate_weights(objects, kind, j)
        if (nrow(ww)) {ww[, constant := constant]; ensemble_weights[[length(ensemble_weights) + 1L]] <- ww}
      }
    }
    if (!all_correct) {
      wrong_ratio[[j]] <- fit_ratios_without_history(task, j, g, shared_splits, nested)
      wrong_ratio[[j]]$population <- predict_ratios(wrong_ratio[[j]]$ratio_fits, pt, g)
    }
    train <- task$folds[[j]]$training_set
    for (s in seq_len(g$tau)) {
      H <- task$vars$history("A", s); A <- task$vars$A[[s]]
      B <- cell_key(task$natural[train, c(H, A), drop = FALSE])
      V <- cell_key(cbind(task$natural[train, H, drop = FALSE],
        d[train, args$health[[s]], drop = FALSE]))
      M <- d[[args$measurement[s]]][train]
      for (kind in c("B", "V")) {
        keys <- if (kind == "B") B else V
        counts <- data.table(cell = keys, M = M)[, .(training = .N, measured = sum(M)), by = cell]
        counts[, `:=`(fold = j, horizon = s, kind = kind)]
        occupancy[[length(occupancy) + 1L]] <- counts
      }
    }
  }
  results <- diagnostics <- functions <- errors <- eif_components <- penalty_failures <- list()
  for (beta_correct in if (all_correct) TRUE else c(TRUE, FALSE))
    for (ratio_correct in if (all_correct) TRUE else c(TRUE, FALSE))
    for (m_correct in if (all_correct) TRUE else c(TRUE, FALSE)) for (method in c("sdr", "tmle")) {
      warnings <- character(); b <- if (beta_correct) 1L else 2L
      callback <- function(training_data, validation_data, inner_folds, task, fold) {
        got <- fitted[[fold]][[b]]
        got$lambda_validation <- fitted[[fold]][[1L]]$lambda_validation
        if (!ratio_correct) {
          got$ratios_training <- wrong_ratio[[fold]]$ratios_training
          got$ratios_validation <- wrong_ratio[[fold]]$ratios_validation
        }
        got
      }
      # A scenario-specific stream keeps subsequent scenarios reproducible
      # if a numerical failure occurs in an earlier one.
      scenario_seed <- seed + 100000L * (as.integer(!beta_correct) * 8L +
        as.integer(!ratio_correct) * 4L + as.integer(!m_correct) * 2L + as.integer(method == "tmle"))
      set.seed(scenario_seed)
      if (is.function(nuisance_callback)) cat(format(Sys.time()), "Fitting sequential estimator", method, "\n")
      value <- tryCatch(withCallingHandlers({
        scenario_args <- args
        scenario_args$bridge_library <- cm_library(!beta_correct, TRUE)
        scenario_args$control$.ratio_history <- ratio_correct
        fit <- do.call(lmtp_bridge, c(list(data = d, estimator = method,
          folds = task$folds, learner_folds = task$learner_folds, nuisance_fit = callback,
          learners_outcome = if (m_correct) sl_library else "SL.time_mean", learners_trt = trt_library), scenario_args))
        q_pop <- vector("list", length(task$folds))
        q_valid <- vector("list", length(task$folds))
        for (j in seq_along(task$folds)) {
          q_pop[[j]] <- predict_regressions(fit$fits_m[[j]], pt, g, method, fit$targeting[[j]])
          valid <- task$folds[[j]]$validation_set
          vt <- study_task(d[valid, ], g)$task
          q_valid[[j]] <- predict_regressions(fit$fits_m[[j]], vt, g, method, fit$targeting[[j]])
          tt <- fit_tuning(fit$fits_m[[j]], "regression", j)
          if (nrow(tt)) {
            tt[, `:=`(beta_correct = beta_correct, ratio_correct = ratio_correct,
                       m_correct = m_correct, estimator = method)]
            tuning[[length(tuning) + 1L]] <- tt
          }
          ww <- fit_candidate_weights(fit$fits_m[[j]], "regression", j)
          if (nrow(ww)) {
            ww[, `:=`(beta_correct = beta_correct, ratio_correct = ratio_correct,
                       m_correct = m_correct)]
            set(ww, j = "estimator", value = method)
            ensemble_weights[[length(ensemble_weights) + 1L]] <- ww
          }
        }
        for (lambda_correct in if (all_correct) TRUE else c(TRUE, FALSE)) {
          l <- if (lambda_correct) 1L else 2L
          scores <- matrix(NA_real_, n, g$tau); seq_scores <- scores
          scenario <- paste0("beta", as.integer(beta_correct), "_lambda", as.integer(lambda_correct),
            "_ratio", as.integer(ratio_correct), "_m", as.integer(m_correct))
          for (j in seq_along(task$folds)) {
            v <- task$folds[[j]]$validation_set
            nu <- c(list(beta = fitted[[j]][[b]]$beta_validation,
              lambda = fitted[[j]][[l]]$lambda_validation,
              ratios = if (ratio_correct) fitted[[j]][[1L]]$ratios_validation else wrong_ratio[[j]]$ratios_validation), q_valid[[j]])
            direct <- raw_scores(d[v, ], g, nu)
            scores[v, ] <- direct$score; seq_scores[v, ] <- direct$sequential
            nu_pop <- c(list(beta = population[[j]][[b]]$beta,
              lambda = population[[j]][[l]]$lambda,
              ratios = if (ratio_correct) population[[j]][[1L]]$ratios else wrong_ratio[[j]]$population), q_pop[[j]])
            dg <- population_diagnostics(p, pt, g, nu_pop, truth_nu, truth)
            dg[, `:=`(fold = j, fold_weight = length(v) / n, scenario = scenario, estimator = method)]
            diagnostics[[length(diagnostics) + 1L]] <- dg
            # Values on the entire finite support fully represent the fitted
            # nuisance functions needed to reproduce these population checks.
            functions[[paste(method, scenario, j, sep = "/")]] <- nu_pop
          }
          if (lambda_correct) for (s in seq_len(g$tau)) {
            stopifnot(max(abs(scores[, s] - fit$estimates[[s]]@eif)) < 1e-8,
              abs(mean(seq_scores[, s]) - fit$sequential_estimates[[s]]@x) < 1e-7)
          }
          estimates <- colMeans(scores); ses <- apply(scores, 2, sd) / sqrt(n)
          if (store_scores) {
            validation_fold <- integer(n)
            for (j in seq_along(task$folds)) validation_fold[task$folds[[j]]$validation_set] <- j
            for (s in seq_len(g$tau)) eif_components[[length(eif_components) + 1L]] <- data.table(
              person = seq_len(n), fold = validation_fold, horizon = s, estimator = method,
              scenario = scenario, sequential_score = seq_scores[, s],
              bridge_adjoint = scores[, s] - seq_scores[, s],
              sequential_eif = seq_scores[, s] - estimates[s],
              total_eif = scores[, s] - estimates[s])
          }
          critical <- joint_critical(scores)
          for (s in seq_len(g$tau)) results[[length(results) + 1L]] <- data.table(
            horizon = s, estimator = method, scenario = scenario,
            beta_correct = beta_correct, lambda_correct = lambda_correct,
            ratio_correct = ratio_correct, m_correct = m_correct,
            sufficient_consistency = (beta_correct || lambda_correct) && (ratio_correct || m_correct),
            truth = truth[s], estimate = estimates[s], se = ses[s],
            mean_sequential_eif = mean(seq_scores[, s]) - truth[s],
            mean_bridge_adjoint_eif = mean(scores[, s] - seq_scores[, s]),
            mean_total_eif = estimates[s] - truth[s],
            lower = estimates[s] - qnorm(.975) * ses[s], upper = estimates[s] + qnorm(.975) * ses[s],
            joint_critical = critical, joint_lower = estimates[s] - critical * ses[s],
            joint_upper = estimates[s] + critical * ses[s],
            joint_curve_covered = all(abs(estimates - truth) <= critical * ses),
            warnings = paste(unique(warnings), collapse = " | "), error = NA_character_)
        }
        TRUE
      }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) {
        if (inherits(e, "cmbridge_penalty_boundary_error"))
          penalty_failures[[length(penalty_failures) + 1L]] <<- list(estimator = method,
            beta_correct = beta_correct, ratio_correct = ratio_correct, m_correct = m_correct, condition = e)
        errors[[length(errors) + 1L]] <<- data.table(beta_correct = beta_correct, ratio_correct = ratio_correct,
          m_correct = m_correct, estimator = method, error = conditionMessage(e), warnings = paste(unique(warnings), collapse = " | "))
        FALSE
      })
    }
  candidate_failures <- penalty_trial_failures <- list()
  for (j in seq_along(fitted)) for (b in seq_along(fitted[[j]])) {
    cache <- fitted[[j]][[b]]$response_nuisance_cache$beta
    if (is.null(cache)) next
    for (key in ls(cache, all.names = TRUE)) {
      entry <- get(key, cache, inherits = FALSE); model <- entry$model
      if (is.null(model)) next
      for (failure in model$candidate_failures) candidate_failures[[length(candidate_failures) + 1L]] <-
        data.table(fold = j, bridge_specification = b, fit_key = key, horizon = entry$horizon,
          fitting_n = model$n, candidate = failure$candidate, scope = failure$scope,
          failure_training_n = failure$training_n, error = failure$message)
      paths <- c(list(lapply(model$candidates, `[[`, "penalty_cv")), model$fold_penalty_cv)
      for (path in paths) for (candidate in names(path)) {
        cv <- path[[candidate]]
        if (!is.null(cv) && "failed" %in% names(cv) && any(cv$failed))
          penalty_trial_failures[[length(penalty_trial_failures) + 1L]] <- cbind(
            as.data.table(cv)[failed == TRUE], data.table(fold = j, bridge_specification = b,
              fit_key = key, horizon = entry$horizon, fitting_n = model$n, candidate = candidate))
      }
    }
  }
  add_design <- function(x) {
    if (!nrow(x)) return(x)
    x[, `:=`(mechanism = mechanism, n = n, replicate = replicate, seed = seed)]
    x
  }
  list(results = add_design(rbindlist(results, fill = TRUE)),
    diagnostics = add_design(rbindlist(diagnostics, fill = TRUE)),
    errors = add_design(rbindlist(errors, fill = TRUE)), penalty_failures = penalty_failures, tuning = rbindlist(tuning, fill = TRUE),
    candidate_failures = add_design(rbindlist(candidate_failures, fill = TRUE)),
    penalty_trial_failures = add_design(rbindlist(penalty_trial_failures, fill = TRUE)),
    ensemble_weights = rbindlist(ensemble_weights, fill = TRUE),
    occupancy = rbindlist(occupancy), population = p, fitted_functions = functions,
    eif_components = add_design(rbindlist(eif_components, fill = TRUE)),
    data = if (store_scores) d else NULL,
    folds = task$folds, learner_folds = task$learner_folds,
    design = list(mechanism = mechanism, n = n, replicate = replicate, seed = seed,
                  all_correct = all_correct, dose_max = g$dose_max, learner_groups = learner_groups,
                  balance_measurement = balance_measurement, visits = visits))
}
