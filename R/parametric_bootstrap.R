parametric_bootstrap <- function(object,
                                      nsim = 200,
                                      alpha = 0.05,
                                      which.pars = "B",
                                      label.method = c("tau", "pi", "hard"),
                                      seed = NULL,
                                      verbose = TRUE) {

  if (!inherits(object, "csam")) {
    stop("object must be of class 'csam'")
  }

  label.method <- match.arg(label.method)

  if (!is.null(seed)) set.seed(seed)

  # Extract model components
  call.list <- as.list(object$call)
  call.list$X <- X <- object$X
  call.list$family <- family <- object$family
  call.list$start.pars <- eval(call.list$start.pars,envir = object$call.env)
  call.list$trace = FALSE
  call.list$verbose = FALSE

  tau_hat <- object$tau

  n <- nrow(object$Y)
  s <- ncol(object$Y)
  p <- ncol(X)
  g <- object$g
  d <- object$d

  pars <- object$par.list
  B_hat <- object$B

  # --- helper: generate permutations ---
  all_perms <- function(g) {
    if (g == 1) return(list(1))
    perms <- gtools::permutations(g, g)
    split(perms, seq(nrow(perms)))
  }

  perms <- all_perms(g)

  # --- helper: align B ---
  align_B <- function(B_boot, B_ref) {
    best_dist <- Inf
    best_B <- B_boot

    for (perm in perms) {
      Bp <- B_boot[perm, , drop = FALSE]
      dist <- sum((Bp - B_ref)^2)
      if (dist < best_dist) {
        best_dist <- dist
        best_B <- Bp
      }
    }
    best_B
  }

  # --- simulate dataset ---
  simulate_from_fit <- function() {

    # sample latent labels
    z <- integer(s)
    for (j in 1:s) {
      if (label.method == "pi") {
        z[j] <- sample(1:g, 1, prob = pars$pi)
      } else if (label.method == "tau") {
        z[j] <- sample(1:g, 1, prob = tau_hat[j, ])
      } else if (label.method == "hard") {
        z[j] <- which.max(tau_hat[j, ])
      }
    }

    eta_g <- tcrossprod(X, pars$B)
    eta <- matrix(NA, n, s)

    for (j in 1:s) {
      eta[, j] <- pars$beta0[j] +
        eta_g[, z[j]] +
        pars$U %*% pars$Lambda[j, ]
    }

    mu <- family$linkinv(eta)

    # response simulation
    Ysim <- matrix(NA, n, s)

    if (family$family == "poisson") {
      Ysim <- matrix(rpois(n * s, mu), n, s)
    } else if (family$family == "binomial") {
      Ysim <- matrix(rbinom(n * s, 1, mu), n, s)
    } else if (family$family == "gaussian") {
      for (j in 1:s) {
        Ysim[, j] <- rnorm(n, mu[, j], sqrt(pars$phi[j]))
      }
    } else if (family$family == "Gamma") {
      for (j in 1:s) {
        shape <- 1 / pars$phi[j]
        scale <- mu[, j] / shape
        Ysim[, j] <- rgamma(n, shape = shape, scale = scale)
      }
    } else {
      stop("Family not implemented")
    }

    Ysim
  }

  # storage
  B_boot <- array(NA, dim = c(nsim, g, p))

  # bootstrap loop
  for (b in 1:nsim) {

    if (verbose && b %% 10 == 0) {
      cat("Bootstrap iteration:", b, "\n")
    }

    call.list$Y <- simulate_from_fit()

    fit_b <- try(
      do.call(csam, call.list[-1]),
      silent = TRUE
    )

    if (inherits(fit_b, "try-error")) next

    B_est <- fit_b$B

    # --- ALIGNMENT STEP ---
    B_est <- align_B(B_est, B_hat)

    B_boot[b, , ] <- B_est
  }

  # remove failed fits or fits where estimates are the result of a spurious fit
  valid <- apply(B_boot, 1, function(x) !any(is.na(x) & all(x < 1e6)))
  B_boot <- B_boot[valid, , , drop = FALSE]

  # CI
  lower <- alpha / 2
  upper <- 1 - alpha / 2

  B_ci_lower <- apply(B_boot, c(2, 3), quantile, probs = lower)
  B_ci_upper <- apply(B_boot, c(2, 3), quantile, probs = upper)

  list(
    B_hat = B_hat,
    B_boot = B_boot,
    CI = list(
      lower = B_ci_lower,
      upper = B_ci_upper,
      level = 1 - alpha
    ),
    label.method = label.method,
    nsim = nsim,
    n_valid = sum(valid)
  )
}
