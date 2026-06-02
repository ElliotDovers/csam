##### SAM EM algorithm components #####

#' E-step posterior probabilities for vanilla SAM (internal)
#'
#' @keywords internal
estep0_post_probs <- function(Y, X, par.list, family,
                              trunc.interval = FALSE) {

  beta0 <- par.list$beta0
  B     <- par.list$B
  pi    <- par.list$pi
  phi   <- par.list$phi   # kept for API compatibility

  n <- nrow(Y)
  s <- ncol(Y)
  g <- nrow(B)

  ## ---- Shared linear predictor ----
  XB <- X %*% t(B)   # n × g

  loglik_mat <- matrix(0, s, g)

  ## ---- E-step: species-level loop (parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    eta_base <- beta0[j]

    for (k in seq_len(g)) {

      eta <- eta_base + XB[, k]
      mu  <- family$linkinv(eta)

      # dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      # loglik_mat[j, k] <- -0.5 * sum(dev)
      if (family$family %in% c("binomial", "poisson")) {
        dev = NULL
      } else {
        dev <- family$dev.resids(y = yj, mu = mu, wt = rep(1, n))
      }
      ll_gj <- -0.5 * family$aic(y = yj, n = rep(1, n), mu = mu, wt = rep(1, n), dev = dev)

      ## relative log-likelihood (constants drop out)
      loglik_mat[j, k] <- sum(ll_gj)
    }
  }

  ## ---- Posterior normalisation ----
  logpost <- sweep(loglik_mat, 2, log(pi), "+")
  row_max <- apply(logpost, 1, max)
  tau_mat <- exp(logpost - row_max)
  tau <- tau_mat / rowSums(tau_mat)

  ## ---- Optional truncation (unchanged) ----
  if (trunc.interval) {
    alpha <- (1 - 0.8 * g) / ((0.8 * 2 - g) - 1)
    tau <- apply(tau, 2, function(x) {
      ((2 * alpha * x) - alpha + 1) /
        ((2 * alpha) - (alpha * g) + g)
    })
  }

  tau
}


#' M-step for vanilla SAM: archetype parameters (internal)
#'
#' @keywords internal
mstep0_arch_pars <- function(Y, X, par.list, tau,
                             family, maxit = 1, use.starts = TRUE) {

  beta0 = par.list$beta0; B = par.list$B; pi = par.list$pi; phi = par.list$phi
  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X)
  B_new <- B

  ## ---- Pre-build shared stacks ONCE ----
  y_stack  <- as.vector(Y)                       # length n*s
  X_stack  <- X[rep(seq_len(n), times = s), ]    # n*s × p

  ns <- n * s
  offset_stack  <- numeric(ns)
  weights_stack <- numeric(ns)

  for (k in seq_len(g)) {

    ## ---- Fill species-wise blocks ----
    for (j in seq_len(s)) {
      idx <- ((j - 1) * n + 1):(j * n)
      offset_stack[idx]  <- beta0[j]
      weights_stack[idx] <- tau[j, k]
    }

    # # note suppressWarnings removes non-convergence as we are typically using only a single iteration
    # fit <- suppressWarnings(stats::glm.fit(
    #   x = X_stack,
    #   y = y_stack,
    #   weights = weights_stack,
    #   offset = offset_stack,
    #   family = family,
    #   control = list(maxit = maxit)
    # ))
    fit <- stats::glm.fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      start = if (use.starts) {B_new[k, ]} else {NULL},
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    )

    if (exists("fit")) {
      B_new[k, ] <- fit$coefficients
    }

  }

  B_new
}

#' M-step for vanilla SAM: species parameters (internal)
#'
#' @keywords internal
mstep0_species_pars <- function(Y, X, par.list, tau,
                                family, maxit = 1, use.starts = TRUE) {

  beta0 = par.list$beta0; B = par.list$B; pi = par.list$pi; phi = par.list$phi
  n <- nrow(Y); s <- ncol(Y); g <- nrow(B)

  beta0_new <- beta0
  phi_new <- phi

  ## ---- Shared linear predictor (fixed in this CM-step) ----
  XB <- X %*% t(B)   # n × g

  ng <- n * g

  ## ---- Allocate stacks ONCE ----
  y_stack       <- numeric(ng)
  offset_stack  <- numeric(ng)
  weights_stack <- numeric(ng)
  X_stack       <- matrix(1, nrow = ng, ncol = 1)  # intercept only

  for (j in seq_len(s)) {

    yj <- Y[, j]

    for (k in seq_len(g)) {
      idx <- ((k - 1) * n + 1):(k * n)
      y_stack[idx]       <- yj
      offset_stack[idx]  <- XB[, k]
      weights_stack[idx] <- tau[j, k]
    }

    # # note suppressWarnings removes non-convergence as we are typically using only a single iteration
    # fit <- suppressWarnings(stats::glm.fit(
    #   x = X_stack,
    #   y = y_stack,
    #   weights = weights_stack,
    #   offset = offset_stack,
    #   family = family,
    #   control = list(maxit = maxit)
    # ))
    fit <- stats::glm.fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      start = if (use.starts) {beta0_new[j]} else {NULL},
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    )

    if (exists("fit")) {
      coef_j <- fit$coefficients
      beta0_new[j] <- coef_j[1]
    }

    if (family$family == "gaussian") {
      mu <- fit$mu
      phi_new[j] <- sum(weights_stack * (y_stack - mu)^2) / sum(weights_stack)
    }
  }

  list(beta0 = beta0_new, phi = phi_new)
}
