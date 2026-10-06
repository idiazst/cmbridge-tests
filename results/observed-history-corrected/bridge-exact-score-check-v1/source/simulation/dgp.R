# Finite-support mechanisms for the observed-history estimator in Section 3.
# Names identify components of C and A; the report uses only paper notation.
# A no-visit treatment is the known deterministic zero prescription.

make_mechanism <- function(mechanism = c("binary_point", "binary_longitudinal", "discrete_dose"),
                           visits = c("scheduled", "intermittent", "missing_both"), strength = 3.6, dose_max = 2L) {
  mechanism <- match.arg(mechanism)
  visits <- match.arg(visits)
  if (length(dose_max) != 1L || is.na(dose_max) || !dose_max %in% c(2L, 3L))
    stop("dose_max must be 2 or 3.")
  tau <- if (mechanism == "binary_point") 1L else 2L
  dose <- mechanism == "discrete_dose"
  missing_both <- visits == "missing_both"
  C <- expand.grid(covariate = 0:1, outcome = 0:1)
  A <- expand.grid(policy_component = if (dose) 0:dose_max else 0:1, other_component = 0:1)
  L <- rbind(data.frame(R = 0, covariate = 0, outcome = 0), cbind(R = 1, C))
  g <- gd <- array(0, c(2, tau, nrow(L), nrow(A)))
  health <- array(0, c(2, tau, nrow(L), nrow(A), nrow(C)))
  measurement <- array(0, c(2, tau, nrow(L), nrow(C)))
  transition <- array(0, c(2, tau, nrow(L), nrow(A), nrow(L)))
  for (baseline in 0:1) for (t in seq_len(tau)) for (state in seq_len(nrow(L))) {
    current <- L[state, ]
    other_p <- plogis(-.15 + .3 * current$outcome - .2 * baseline)
    if (dose) {
      log_weights <- c(0, .15 + .3 * current$covariate - .2 * current$outcome,
                      -.1 + .4 * current$covariate - .3 * current$outcome + .1 * baseline)
      if (dose_max == 3L) log_weights <- c(log_weights,
        -.35 + .5 * current$covariate - .4 * current$outcome + .2 * baseline)
      first_p <- exp(log_weights) / sum(exp(log_weights))
    } else {
      probability <- plogis(-.2 + .45 * current$covariate - .4 * current$outcome + .15 * baseline)
      first_p <- c(1 - probability, probability)
    }
    if (current$R == 0) {
      # alpha_t(H_t) = (0,0), a known deterministic prescription rule.
      g[baseline + 1, t, state, 1] <- 1
    } else {
      g[baseline + 1, t, state, ] <- first_p[A$policy_component + 1] *
        dbinom(A$other_component, 1, other_p)
    }
    for (action in seq_len(nrow(A))) {
      changed <- if (current$R == 0) A$policy_component[action] else
        if (dose) min(A$policy_component[action] + 1, dose_max) else 1
      destination <- which(A$policy_component == changed &
                             A$other_component == A$other_component[action])
      gd[baseline + 1, t, state, destination] <- gd[baseline + 1, t, state, destination] +
        g[baseline + 1, t, state, action]
      covariate_p <- plogis(-strength / 2 + strength * A$other_component[action] +
                             .35 * A$policy_component[action] + .3 * current$covariate + .25 * baseline)
      outcome_p <- plogis(-strength / 2 + strength * A$policy_component[action] / (if (dose) 2 else 1) +
                           .35 * C$covariate + .3 * current$outcome + .2 * baseline)
      # At a missed visit A_t is deterministic. The target-specific adjoint
      # then requires the next outcome to be known from H_t whenever the
      # preceding policy weight is positive. Retain the baseline outcome.
      if (missing_both && current$R == 0) outcome_p <- rep(baseline, nrow(C))
      health[baseline + 1, t, state, action, ] <- dbinom(C$covariate, 1, covariate_p) *
        dbinom(C$outcome, 1, outcome_p)
    }
    measurement[baseline + 1, t, state, ] <- if (visits == "scheduled" && t < tau) 1 else
      plogis(.25 - .65 * C$outcome - .45 * C$covariate + .15 * current$outcome + .1 * baseline)
    for (action in seq_len(nrow(A))) {
      transition[baseline + 1, t, state, action, 1] <-
        sum(health[baseline + 1, t, state, action, ] * (1 - measurement[baseline + 1, t, state, ]))
      transition[baseline + 1, t, state, action, 2:5] <-
        health[baseline + 1, t, state, action, ] * measurement[baseline + 1, t, state, ]
    }
  }
  ratio <- array(0, dim(g))
  positive <- g > 0
  ratio[positive] <- gd[positive] / g[positive]
  result <- list(mechanism = mechanism, visits = visits, tau = tau, dose = dose,
                 dose_max = if (dose) dose_max else 1L, C = C, A = A, L = L, g = g, gd = gd, health = health,
                 measurement = measurement, transition = transition,
                 ratio = ratio, beta = 1 / measurement,
                 initial_outcomes = if (missing_both) 0:1 else c(0L, 0L))
  if (missing_both) for (baseline in 1:2) for (t in seq_len(tau))
    for (state in which(L$R == 0)) {
      action <- which(g[baseline, t, state, ] > 0)
      stopifnot(length(action) == 1L)
      observed_probability <- sum(health[baseline, t, state, action, ] *
                                    measurement[baseline, t, state, ])
      # One valid, observable bridge solution at this history is constant
      # in next health. No uniqueness is assumed for these histories.
      result$beta[baseline, t, state, ] <- 1 / observed_probability
    }
  result$natural_states <- state_probabilities(result, policy = FALSE)
  result$policy_states <- state_probabilities(result, policy = TRUE)
  result
}

state_probabilities <- function(mechanism, policy = FALSE) {
  value <- array(0, c(2, mechanism$tau + 1L, nrow(mechanism$L)))
  initial_outcomes <- if (is.null(mechanism$initial_outcomes)) c(0L, 0L) else mechanism$initial_outcomes
  for (baseline in 1:2)
    value[baseline, 1, 2L + 2L * initial_outcomes[baseline]] <- .5
  treatment <- if (policy) mechanism$gd else mechanism$g
  for (baseline in 1:2) for (t in seq_len(mechanism$tau)) {
    for (state in seq_len(nrow(mechanism$L))) for (action in seq_len(nrow(mechanism$A))) {
      value[baseline, t + 1, ] <- value[baseline, t + 1, ] + value[baseline, t, state] *
        treatment[baseline, t, state, action] * mechanism$transition[baseline, t, state, action, ]
    }
  }
  value
}

exact_truth <- function(mechanism) {
  answer <- numeric(mechanism$tau)
  for (baseline in 1:2) for (t in seq_len(mechanism$tau)) {
    for (state in seq_len(nrow(mechanism$L))) for (action in seq_len(nrow(mechanism$A))) {
      answer[t] <- answer[t] + mechanism$policy_states[baseline, t, state] *
        mechanism$gd[baseline, t, state, action] *
        sum(mechanism$health[baseline, t, state, action, ] * mechanism$C$outcome)
    }
  }
  answer
}

backward_truth <- function(mechanism) {
  answer <- numeric(mechanism$tau)
  for (endpoint in seq_len(mechanism$tau)) {
    value <- matrix(0, 2, nrow(mechanism$L))
    for (t in endpoint:1) {
      next_value <- value
      for (baseline in 1:2) for (state in seq_len(nrow(mechanism$L))) {
        conditional <- vapply(seq_len(nrow(mechanism$A)), function(action) {
          if (t == endpoint) sum(mechanism$health[baseline, t, state, action, ] * mechanism$C$outcome) else
            sum(mechanism$transition[baseline, t, state, action, ] * next_value[baseline, ])
        }, numeric(1))
        value[baseline, state] <- sum(mechanism$gd[baseline, t, state, ] * conditional)
      }
    }
    answer[endpoint] <- sum(mechanism$natural_states[, 1, ] * value)
  }
  answer
}

audit_mechanism <- function(mechanism, tolerance = 1e-10) {
  rows <- list()
  for (baseline in 1:2) for (t in seq_len(mechanism$tau)) {
    for (state in seq_len(nrow(mechanism$L))) {
      natural <- mechanism$natural_states[baseline, t, state]
      policy <- mechanism$policy_states[baseline, t, state]
      if (natural + policy <= tolerance) next
      probability <- mechanism$g[baseline, t, state, ]
      active <- which(probability > 0)
      K <- mechanism$health[baseline, t, state, , , drop = TRUE]
      joint <- K * probability
      marginal <- colSums(joint)
      C_active <- which(marginal > 0)
      posterior <- t(joint[, C_active, drop = FALSE]) / marginal[C_active]
      loading <- mechanism$C$outcome[C_active] *
        as.vector(posterior %*% mechanism$ratio[baseline, t, state, ])
      decomposition <- svd(posterior)
      keep <- which(decomposition$d > tolerance * max(decomposition$d))
      adjoint <- if (length(keep)) as.vector(decomposition$v[, keep, drop = FALSE] %*%
        (as.vector(crossprod(decomposition$u[, keep, drop = FALSE], loading)) / decomposition$d[keep])) else
        rep(0, nrow(mechanism$A))
      dual_error <- max(abs(as.vector(posterior %*% adjoint) - loading))
      operator <- sweep(K[active, C_active, drop = FALSE], 2,
                        mechanism$measurement[baseline, t, state, C_active], `*`)
      spectrum <- svd(operator)$d
      bridge_error <- max(abs(as.vector(operator %*% mechanism$beta[baseline, t, state, C_active]) - 1))
      no_visit <- mechanism$L$R[state] == 0
      rows[[length(rows) + 1L]] <- data.frame(
        mechanism = mechanism$mechanism, visits = mechanism$visits, time = t,
        C1_baseline = baseline - 1L, recorded_state = state,
        natural_probability = natural, policy_probability = policy,
        no_visit = no_visit, treatment_values = length(active),
        assumption_1 = !no_visit || (length(active) == 1L && active == 1L),
        assumption_3 = all(mechanism$gd[baseline, t, state, probability == 0] == 0) &&
          !(policy > tolerance && natural <= tolerance),
        assumption_5 = min(mechanism$measurement[baseline, t, state, C_active]) > 0,
        assumption_7 = dual_error < tolerance,
        measurement_min = min(mechanism$measurement[baseline, t, state, C_active]),
        operator_rank = sum(spectrum > tolerance * max(spectrum)),
        minimum_singular_value = min(spectrum),
        bridge_equation_error = bridge_error,
        adjoint_equation_error = dual_error)
    }
  }
  do.call(rbind, rows)
}

require_valid_mechanism <- function(mechanism) {
  stopifnot(all(is.finite(mechanism$g)), all(is.finite(mechanism$gd)),
            all(is.finite(mechanism$health)), all(is.finite(mechanism$measurement)),
            all(mechanism$g >= 0), all(mechanism$gd >= 0), all(mechanism$health >= 0),
            all(mechanism$measurement > 0 & mechanism$measurement <= 1),
            max(abs(apply(mechanism$g, 1:3, sum) - 1)) < 1e-12,
            max(abs(apply(mechanism$gd, 1:3, sum) - 1)) < 1e-12,
            max(abs(apply(mechanism$health, 1:4, sum) - 1)) < 1e-12,
            max(abs(apply(mechanism$transition, 1:4, sum) - 1)) < 1e-12)
  audit <- audit_mechanism(mechanism)
  required <- c("assumption_1", "assumption_3", "assumption_5", "assumption_7")
  if (!all(as.matrix(audit[, required]))) {
    failed <- required[!vapply(audit[, required], all, logical(1))]
    stop("Mechanism fails required checks: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  invisible(audit)
}

draw_data <- function(n, mechanism) {
  # A failed population audit prevents even a pilot sample being generated.
  require_valid_mechanism(mechanism)
  result <- data.frame(C1_baseline = rbinom(n, 1, .5), C1_covariate = 0, Y1 = 0, R1 = 1)
  if (mechanism$visits == "missing_both") result$Y1 <- result$C1_baseline
  categories <- function(probability) {
    uniform <- runif(nrow(probability))
    value <- rep(1L, nrow(probability))
    cumulative <- probability[, 1]
    for (j in 2:ncol(probability)) {
      value[uniform > cumulative] <- j
      cumulative <- cumulative + probability[, j]
    }
    value
  }
  for (t in seq_len(mechanism$tau)) {
    state <- ifelse(result[[paste0("R", t)]] == 0, 1L,
      2L + ifelse(is.na(result[[paste0("C", t, "_covariate")]]), 0, result[[paste0("C", t, "_covariate")]]) +
        2L * ifelse(is.na(result[[paste0("Y", t)]]), 0, result[[paste0("Y", t)]]))
    probability <- sapply(seq_len(nrow(mechanism$A)), function(a)
      mechanism$g[cbind(result$C1_baseline + 1L, t, state, a)])
    action <- categories(probability)
    result[[paste0("A", t, "_policy")]] <- mechanism$A$policy_component[action]
    result[[paste0("A", t, "_other")]] <- mechanism$A$other_component[action]
    probability <- sapply(seq_len(nrow(mechanism$C)), function(c)
      mechanism$health[cbind(result$C1_baseline + 1L, t, state, action, c)])
    complete <- categories(probability)
    recorded <- rbinom(n, 1,
      mechanism$measurement[cbind(result$C1_baseline + 1L, t, state, complete)])
    result[[paste0("R", t + 1L)]] <- recorded
    result[[paste0("C", t + 1L, "_covariate")]] <- ifelse(recorded == 1, mechanism$C$covariate[complete], NA_real_)
    result[[paste0("Y", t + 1L)]] <- ifelse(recorded == 1, mechanism$C$outcome[complete], NA_real_)
  }
  result
}
