#' Penalized GLM fit (internal)
#'
#' @keywords internal
penalized_glm_fit <- function(X, y, weights = NULL, offset = NULL,
                              penalty_weights = NULL, lambda = 0,
                              family = gaussian(), start = NULL,
                              maxit = 50, tol = 1e-6, verbose = FALSE) {

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset)) offset <- rep(0, n)
  if (is.null(penalty_weights)) penalty_weights <- rep(0, p)

  P <- diag(penalty_weights, p, p)

  linkinv  <- family$linkinv
  mu_eta   <- family$mu.eta
  variance <- family$variance

  if (is.null(start)) {
    beta <- rep(0, p)
  } else {
    if (length(start) == p) {
      beta <- start
    } else {
      stop("incorrect length of starting values")
    }
  }

  eta <- as.vector(X %*% beta + offset)
  mu  <- linkinv(eta)

  for (iter in 1:maxit) {

    mu_eta_val <- mu_eta(eta)
    var_mu <- variance(mu)

    W <- weights * (mu_eta_val^2 / var_mu)
    z <- (eta - offset) + (y - mu) / mu_eta_val

    WX <- sqrt(W) * X
    wz <- sqrt(W) * z

    A <- rbind(WX, sqrt(lambda) * P)
    b <- c(wz, rep(0, p))

    qrA <- qr(A)
    beta_new <- qr.coef(qrA, b)

    eta_new <- as.vector(X %*% beta_new + offset)
    mu_new  <- linkinv(eta_new)

    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new; eta <- eta_new; mu <- mu_new
      break
    }

    beta <- beta_new
    eta <- eta_new
    mu  <- mu_new
  }

  list(coefficients = beta, eta = eta, mu = mu)
}

#' E-step posterior probabilities (internal)
#'
#' @keywords internal
estep_post_probs <- function(Y, X, beta0, B, U, Lambda, phi, pi, family,
                             trunc.interval = FALSE) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B)

  loglik_mat <- matrix(NA, s, g)

  for (j in 1:s) {
    yj <- Y[, j]
    lambda_j <- Lambda[j, , drop = FALSE]
    offset_common <- as.vector(U %*% t(lambda_j))

    for (k in 1:g) {
      eta <- beta0[j] + X %*% B[k, ] + offset_common
      mu  <- family$linkinv(eta)

      dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      ll  <- -0.5 * family$aic(yj, rep(1, n), mu, rep(1, n), dev)

      loglik_mat[j, k] <- sum(ll)
    }
  }

  logpost <- sweep(loglik_mat, 2, log(pi), "+")
  tau <- exp(logpost - apply(logpost, 1, max))
  tau <- tau / rowSums(tau)

  if (trunc.interval) {
    alpha <- (1 - 0.8 * g) / ((0.8 * 2 - g) - 1)
    tau <- apply(tau, 2, function(x) {
      ((2 * alpha * x) - alpha + 1) / ((2 * alpha) - (alpha * g) + g)
    })
  }

  tau
}

#' M-step: archetype parameters (internal)
#'
#' @keywords internal
mstep_arch_pars <- function(Y, X, beta0, B, U, Lambda, phi, pi, tau,
                            family, maxit = 1) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X)
  B_new <- B

  for (k in 1:g) {

    y_stack <- as.vector(Y)
    X_stack <- do.call(rbind, replicate(s, X, simplify = FALSE))

    offset_stack <- rep(0, n * s)
    weights_stack <- rep(0, n * s)

    for (j in 1:s) {
      idx <- ((j - 1) * n + 1):(j * n)
      offset_stack[idx] <- beta0[j] + as.vector(U %*% Lambda[j, ])
      weights_stack[idx] <- tau[j, k]
    }

    # note suppressWarnings removes non-convergence as we are typically using only a single iteration
    fit <- suppressWarnings(stats::glm.fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    ))

    B_new[k, ] <- fit$coefficients
  }

  B_new
}

#' M-step: species parameters (internal)
#'
#' @keywords internal
mstep_species_pars <- function(Y, X, B, U, beta0, Lambda, phi, pi, tau,
                               family, psi2 = 0, maxit = 1) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); d <- ncol(U)

  beta0_new <- beta0
  Lambda_new <- Lambda
  phi_new <- phi

  for (j in 1:s) {

    y_stack <- rep(NA, n * g)
    X_stack <- matrix(0, nrow = n * g, ncol = 1 + d)
    weights_stack <- rep(0, n * g)
    offset_stack <- rep(0, n * g)

    for (k in 1:g) {
      idx <- ((k - 1) * n + 1):(k * n)
      y_stack[idx] <- Y[, j]
      X_stack[idx, 1] <- 1
      X_stack[idx, -1] <- U
      offset_stack[idx] <- as.vector(X %*% B[k, ])
      weights_stack[idx] <- tau[j, k]
    }

    penalty_weights <- c(0, rep(1, d))

    fit <- penalized_glm_fit(
      X = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      penalty_weights = penalty_weights,
      lambda = psi2,
      family = family,
      maxit = maxit,
      start = c(beta0_new[j], Lambda_new[j, ])
    )

    coef_j <- fit$coefficients
    beta0_new[j] <- coef_j[1]
    Lambda_new[j, ] <- coef_j[-1]

    if (family$family == "gaussian") {
      mu <- fit$mu
      phi_new[j] <- sum(weights_stack * (y_stack - mu)^2) / sum(weights_stack)
    }
  }

  list(beta0 = beta0_new, Lambda = Lambda_new, phi = phi_new)
}

#' M-step: site scores (internal)
#'
#' @keywords internal
mstep_site_scores <- function(Y, X, B, beta0, Lambda, U, phi, pi, tau,
                              family, psi1 = 0, maxit = 1) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); d <- ncol(Lambda)
  U_new <- U

  for (i in 1:n) {

    y_stack <- rep(NA, s * g)
    X_stack <- matrix(0, nrow = s * g, ncol = d)
    weights_stack <- rep(0, s * g)
    offset_stack <- rep(0, s * g)

    row <- 1
    for (j in 1:s) {
      for (k in 1:g) {
        y_stack[row] <- Y[i, j]
        X_stack[row, ] <- Lambda[j, ]
        offset_stack[row] <- beta0[j] + as.numeric(X[i, ] %*% B[k, ])
        weights_stack[row] <- tau[j, k]
        row <- row + 1
      }
    }

    penalty_weights <- rep(1, d)

    fit <- penalized_glm_fit(
      X = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      penalty_weights = penalty_weights,
      lambda = psi1,
      family = family,
      maxit = maxit,
      start = U_new[i, ]
    )

    U_new[i, ] <- fit$coefficients
  }

  U_new
}

#' E-step posterior probabilities for vanilla SAM (internal)
#'
#' @keywords internal
estep0_post_probs <- function(Y, X, beta0, B, phi, pi, family,
                             trunc.interval = FALSE) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B)

  loglik_mat <- matrix(NA, s, g)

  for (j in 1:s) {
    yj <- Y[, j]

    for (k in 1:g) {
      eta <- beta0[j] + X %*% B[k, ]
      mu  <- family$linkinv(eta)

      dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      ll  <- -0.5 * family$aic(yj, rep(1, n), mu, rep(1, n), dev)

      loglik_mat[j, k] <- sum(ll)
    }
  }

  logpost <- sweep(loglik_mat, 2, log(pi), "+")
  tau <- exp(logpost - apply(logpost, 1, max))
  tau <- tau / rowSums(tau)

  if (trunc.interval) {
    alpha <- (1 - 0.8 * g) / ((0.8 * 2 - g) - 1)
    tau <- apply(tau, 2, function(x) {
      ((2 * alpha * x) - alpha + 1) / ((2 * alpha) - (alpha * g) + g)
    })
  }

  tau
}

#' M-step for vanilla SAM: archetype parameters (internal)
#'
#' @keywords internal
mstep0_arch_pars <- function(Y, X, beta0, B, phi, pi, tau,
                            family, maxit = 1) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X)
  B_new <- B

  for (k in 1:g) {

    y_stack <- as.vector(Y)
    X_stack <- do.call(rbind, replicate(s, X, simplify = FALSE))

    offset_stack <- rep(0, n * s)
    weights_stack <- rep(0, n * s)

    for (j in 1:s) {
      idx <- ((j - 1) * n + 1):(j * n)
      offset_stack[idx] <- beta0[j]
      weights_stack[idx] <- tau[j, k]
    }

    # note suppressWarnings removes non-convergence as we are typically using only a single iteration
    fit <- suppressWarnings(stats::glm.fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    ))

    B_new[k, ] <- fit$coefficients
  }

  B_new
}

#' M-step for vanilla SAM: species parameters (internal)
#'
#' @keywords internal
mstep0_species_pars <- function(Y, X, B, beta0, phi, pi, tau,
                               family, maxit = 1) {

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B)

  beta0_new <- beta0
  phi_new <- phi

  for (j in 1:s) {

    y_stack <- rep(NA, n * g)
    X_stack <- matrix(0, nrow = n * g, ncol = 1)
    weights_stack <- rep(0, n * g)
    offset_stack <- rep(0, n * g)

    for (k in 1:g) {
      idx <- ((k - 1) * n + 1):(k * n)
      y_stack[idx] <- Y[, j]
      X_stack[idx, 1] <- 1
      offset_stack[idx] <- as.vector(X %*% B[k, ])
      weights_stack[idx] <- tau[j, k]
    }

    # note suppressWarnings removes non-convergence as we are typically using only a single iteration
    fit <- suppressWarnings(stats::glm.fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    ))

    coef_j <- fit$coefficients
    beta0_new[j] <- coef_j[1]

    if (family$family == "gaussian") {
      mu <- fit$mu
      phi_new[j] <- sum(weights_stack * (y_stack - mu)^2) / sum(weights_stack)
    }
  }

  list(beta0 = beta0_new, phi = phi_new)
}


#' Calculate the complete-data expected score function (internal)
#'
#' @keywords internal
score_complete <- function(object) {

  Y <- object$Y; X <- object$X; family <- object$family
  beta0 <- object$beta0; B <- object$B; pi <- object$pi
  U <- object$U; Lambda <- object$Lambda; phi <- object$phi
  psi1 <- object$psi1; psi2 <- object$psi2

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)
  g <- nrow(B); d <- ncol(Lambda)
  is_gaussian <- identical(family$family, "gaussian")

  p_theta <- s + g*p + (g-1) + n*d + s*d + if (is_gaussian) s else 0

  S <- matrix(0, s*g, p_theta)
  row_id <- 1

  for (j in 1:s) {
    for (k in 1:g) {

      lambda_j <- Lambda[j, , drop=FALSE]
      eta <- beta0[j] + X %*% B[k, ] + as.vector(U %*% t(lambda_j))

      mu <- family$linkinv(eta)
      var_mu <- family$variance(mu)
      mu_eta <- family$mu.eta(eta)
      w <- (1/var_mu) * mu_eta
      resid <- Y[, j] - mu
      wr <- as.numeric(w * resid)   # ensure vector

      score_beta0 <- rep(0, s)
      score_beta0[j] <- sum(wr)

      score_B <- rep(0, g*p)
      B_start <- (k-1)*p + 1
      score_B[B_start:(B_start+p-1)] <- colSums(wr * X)

      if (g > 1) {
        pi_sub <- pi[1:(g-1)]
        score_pi <- -pi_sub
        if (k <= g-1) score_pi[k] <- score_pi[k] + 1
      } else {
        score_pi <- numeric(0)
      }

      score_U <- rep(wr, each = d) * rep(as.numeric(lambda_j), times = n)

      score_Lambda <- rep(0, s*d)
      L_start <- (j-1)*d + 1
      U_weighted <- sweep(U, 1, wr, `*`)
      score_Lambda[L_start:(L_start+d-1)] <- colSums(U_weighted)

      if (is_gaussian) {
        score_phi <- rep(0, s)
        phi_j <- phi[j]
        s_phi <- -0.5 * sum((Y[, j] - mu)^2)/phi_j + 0.5*n
        score_phi[j] <- phi_j * s_phi
      } else {
        score_phi <- numeric(0)
      }

      score_U <- score_U - psi1 * as.vector(t(U))
      score_Lambda <- score_Lambda - psi2 * as.vector(t(Lambda))

      S[row_id, ] <- c(score_beta0, score_B, score_pi, score_U, score_Lambda, score_phi)
      row_id <- row_id + 1
    }
  }

  S
}

#' Calculate the complete-data expected Hessian (internal)
#'
#' @keywords internal
hessian_complete <- function(object) {

  Y <- object$Y; X <- object$X; family <- object$family
  beta0 <- object$beta0; B <- object$B; pi <- object$pi
  U <- object$U; Lambda <- object$Lambda; phi <- object$phi
  psi1 <- object$psi1; psi2 <- object$psi2

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)
  g <- nrow(B); d <- ncol(Lambda)
  is_gaussian <- identical(family$family, "gaussian")

  p_theta <- s + g*p + (g-1) + n*d + s*d + if (is_gaussian) s else 0

  H_list <- vector("list", s*g)
  row_id <- 1

  for (j in 1:s) {
    for (k in 1:g) {

      lambda_j <- Lambda[j, , drop=FALSE]
      eta <- beta0[j] + X %*% B[k, ] + as.vector(U %*% t(lambda_j))

      mu <- family$linkinv(eta)
      var_mu <- family$variance(mu)
      mu_eta <- family$mu.eta(eta)
      w <- (1/var_mu) * mu_eta

      H_jk <- matrix(0, p_theta, p_theta)

      for (i in 1:n) {

        Ji_beta0 <- rep(0, s)
        Ji_beta0[j] <- 1

        Ji_B <- rep(0, g*p)
        B_start <- (k-1)*p + 1
        Ji_B[B_start:(B_start+p-1)] <- X[i, ]

        Ji_pi <- if (g > 1) rep(0, g-1) else numeric(0)

        Ji_U <- rep(0, n*d)
        U_pos <- (i-1)*d + 1
        Ji_U[U_pos:(U_pos+d-1)] <- as.numeric(lambda_j)

        Ji_Lambda <- rep(0, s*d)
        L_start <- (j-1)*d + 1
        Ji_Lambda[L_start:(L_start+d-1)] <- as.numeric(U[i, ])

        Ji_phi <- if (is_gaussian) rep(0, s) else numeric(0)

        Ji <- c(Ji_beta0, Ji_B, Ji_pi, Ji_U, Ji_Lambda, Ji_phi)

        H_jk <- H_jk + w[i] * (Ji %o% Ji)
      }

      if (g > 1) {
        pi_sub <- pi[1:(g-1)]
        Lp <- length(pi_sub)
        S_pi <- -diag(pi_sub, Lp, Lp) + outer(pi_sub, pi_sub)
        pi_start <- s + g*p + 1
        idx <- pi_start:(pi_start+Lp-1)
        H_jk[idx, idx] <- H_jk[idx, idx] - S_pi
      }

      U_start <- s + g*p + max(0, g-1) + 1
      Lambda_start <- U_start + n*d

      if (psi1 > 0) {
        idx <- U_start:(U_start+n*d-1)
        H_jk[idx, idx] <- H_jk[idx, idx] + psi1 * diag(n*d)
      }
      if (psi2 > 0) {
        idx <- Lambda_start:(Lambda_start+s*d-1)
        H_jk[idx, idx] <- H_jk[idx, idx] + psi2 * diag(s*d)
      }

      H_list[[row_id]] <- H_jk
      row_id <- row_id + 1
    }
  }

  H_list
}

#' Initialise starting values for Species Archetype Models
#'
#' Generates starting values for a Species Archetype Model (SAM) by fitting
#' species‑specific GLMs to obtain warm‑start coefficients and clustering the
#' resulting slope estimates into \eqn{g} archetypes using k‑means. The returned
#' list provides initial values for species intercepts, archetype slopes, mixing
#' proportions, and dispersion parameters, suitable for use in ECM optimisation.
#'
#' @param Y A numeric response matrix with \eqn{n} rows (sites) and \eqn{s}
#'   columns (species).
#' @param X A numeric design matrix with \eqn{n} rows and \eqn{p} predictors.
#' @param g Integer specifying the number of archetypes. Defaults to 3.
#' @param family A GLM family object used for the species‑specific warm‑start
#'   models. Defaults to \code{poisson()}.
#'
#' @details
#' The function proceeds in two stages:
#'
#' \strong{1. Species‑specific warm starts.}
#' Each species is fitted with a GLM using \code{glm.fit()}, yielding an
#' intercept and a vector of slopes. These form a matrix of species‑level slope
#' estimates of dimension \eqn{s \times p}.
#'
#' \strong{2. Archetype clustering.}
#' The slope matrix is clustered into \eqn{g} groups using
#' \code{stats::kmeans()}. The cluster centres initialise the archetype‑level
#' slope matrix \eqn{B}, and the empirical cluster frequencies initialise the
#' mixing proportions \eqn{\pi}.
#'
#' The returned list contains:
#' \itemize{
#'   \item \code{beta0}: species‑specific intercepts;
#'   \item \code{B}: a \eqn{g \times p} matrix of archetype slopes;
#'   \item \code{pi}: mixing proportions for the \eqn{g} archetypes;
#'   \item \code{phi}: species‑specific dispersion parameters (initially 1).
#' }
#'
#' An attribute \code{"sp clust"} stores the species‑level cluster assignments
#' from the k‑means step.
#'
#' @return A list containing initial parameter values for a SAM, with an
#'   attribute storing species cluster assignments.
#'
#' @export
#'
#' @importFrom stats kmeans glm.fit
init.sam.pars <- function(Y, X, g = 3, family = poisson()) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise parameter list
  start.pars = list(beta0 = rep(0, s), B = matrix(rep(0, g * p), g, p), pi = rep(1/g, g), phi = rep(1, s))

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013
  beta0.init = vector("numeric", s)
  beta1.init = matrix(rep(0, s * p), s, p)
  for (sp in 1:s) {
    tmp.m = glm.fit(x = cbind(rep(1, nrow(X)), X), y = Y[,sp], family = family)
    beta0.init[sp] = tmp.m$coefficients[1]
    beta1.init[sp, ] = tmp.m$coefficients[2:(p + 1)]
  }
  clust = stats::kmeans(beta1.init, g)
  # update the starting parameters for the slopes and intercepts
  start.pars$beta0 = beta0.init
  start.pars$B = clust$centers
  start.pars$pi = prop.table(table(clust$cluster))

  attr(start.pars, "sp clust") = clust$clust
  start.pars
}

#' Initialise starting values for correlated Species Archetype Models
#'
#' Generates starting values for a correlated Species Archetype Model (CSAM) by fitting
#' species‑specific GLMs to obtain warm‑start coefficients and clustering the
#' resulting slope estimates into \eqn{g} archetypes using k‑means. To obtain starting
#' values for the factor-analytic components, a \eqn{d}-latent factor analysis is performed
#' on the residuals from the above archetypal model (with hard classification). The returned
#' list provides initial values for species intercepts, archetype slopes, mixing
#' proportions, and dispersion parameters, as well as the site-specific factor scores and
#' species-specific factor loadings. Ready for input into the \code{csam} function via the
#' \code{init} argument.
#'
#' @param Y A numeric response matrix with \eqn{n} rows (sites) and \eqn{s}
#'   columns (species).
#' @param X A numeric design matrix with \eqn{n} rows and \eqn{p} predictors.
#' @param g Integer specifying the number of archetypes. Defaults to 3.
#' @param family A GLM family object used for the species‑specific warm‑start
#'   models. Defaults to \code{poisson()}.
#' @param d Integer specifying the number of latent factors. Defaults to 2.
#'
#' @details
#' The function proceeds in three stages:
#'
#' \strong{1. Species‑specific warm starts.}
#' Each species is fitted with a GLM using \code{glm.fit()}, yielding an
#' intercept and a vector of slopes. These form a matrix of species‑level slope
#' estimates of dimension \eqn{s \times p}.
#'
#' \strong{2. Archetype clustering.}
#' The slope matrix is clustered into \eqn{g} groups using
#' \code{stats::kmeans()}. The cluster centres initialise the archetype‑level
#' slope matrix \eqn{B}, and the empirical cluster frequencies initialise the
#' mixing proportions \eqn{\pi}.
#'
#' \strong{3. Residual Factor Analysis.}
#' Randomized quantile residuals from the above \eqn{g} clustering (computed via \code{statmod::qresid}) are run in a \eqn{d} latent factor analysis
#' using \code{gllvm:gllvm()}. Scores and loadings form the matrices of these parameters within the CSAM.
#'
#' The returned list contains:
#' \itemize{
#'   \item \code{beta0}: species‑specific intercepts;
#'   \item \code{B}: a \eqn{g \times p} matrix of archetype slopes;
#'   \item \code{pi}: mixing proportions for the \eqn{g} archetypes;
#'   \item \code{U}: a \eqn{n \times d} matrix of site-specific factor scores;
#'   \item \code{Lambda}: a \eqn{s \times d} matrix of species-specific factor loadings;
#'   \item \code{phi}: species‑specific dispersion parameters (initially 1).
#' }
#'
#' An attribute \code{"sp clust"} stores the species‑level cluster assignments
#' from the k‑means step.
#'
#' @return A list containing initial parameter values for a SAM, with an
#'   attribute storing species cluster assignments.
#'
#' @export
#'
#' @importFrom gllvm gllvm
#' @importFrom statmod qresid
init.fa.pars_gllvm <- function(Y, X, g = 3, family = poisson(), d = 2) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise usual SAM parameters
  start.pars = csam::init.sam.pars(Y, X, g = 3, family = family)

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013
  beta0.init = vector("numeric", s)
  res = matrix(rep(0, n * s), n, s)
  for (sp in 1:s) {
    tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
    tmp.m = glm.fit(x = rep(1, nrow(X)), y = Y[,sp], family = family, offset = tmp.offset)
    res[ , sp] = statmod::qresid(tmp.m)
  }
  # perform factor analysis on the residuals
  tmp.fa = gllvm::gllvm(y = res, X = matrix(rep(1, n), ncol = 1), family = gaussian(), num.lv = d)

  # update the starting parameters for the slopes and intercepts
  start.pars$U = tmp.fa$lvs
  start.pars$Lambda = coef(tmp.fa)$Species.scores

  start.pars[c("beta0", "B", "pi", "U", "Lambda", "phi")]
}

#' Initialise starting values for correlated Species Archetype Models
#'
#' Generates starting values for a correlated Species Archetype Model (CSAM) by fitting
#' species‑specific GLMs to obtain warm‑start coefficients and clustering the
#' resulting slope estimates into \eqn{g} archetypes using k‑means. To obtain starting
#' values for the factor-analytic components, a \eqn{d}-latent factor analysis is performed
#' on the residuals from the above archetypal model (with hard classification). The returned
#' list provides initial values for species intercepts, archetype slopes, mixing
#' proportions, and dispersion parameters, as well as the site-specific factor scores and
#' species-specific factor loadings. Ready for input into the \code{csam} function via the
#' \code{init} argument.
#'
#' @param Y A numeric response matrix with \eqn{n} rows (sites) and \eqn{s}
#'   columns (species).
#' @param X A numeric design matrix with \eqn{n} rows and \eqn{p} predictors.
#' @param g Integer specifying the number of archetypes. Defaults to 3.
#' @param family A GLM family object used for the species‑specific warm‑start
#'   models. Defaults to \code{poisson()}.
#' @param d Integer specifying the number of latent factors. Defaults to 2.
#'
#' @details
#' The function proceeds in three stages:
#'
#' \strong{1. Species‑specific warm starts.}
#' Each species is fitted with a GLM using \code{glm.fit()}, yielding an
#' intercept and a vector of slopes. These form a matrix of species‑level slope
#' estimates of dimension \eqn{s \times p}.
#'
#' \strong{2. Archetype clustering.}
#' The slope matrix is clustered into \eqn{g} groups using
#' \code{stats::kmeans()}. The cluster centres initialise the archetype‑level
#' slope matrix \eqn{B}, and the empirical cluster frequencies initialise the
#' mixing proportions \eqn{\pi}.
#'
#' \strong{3. Residual Factor Analysis.}
#' Randomized quantile residuals from the above \eqn{g} clustering (computed via \code{residuals.glmmTMB(x, type = "dunn-smyth")}) are run in a \eqn{d}-latent factor analysis
#' using \code{glmmTMB}'s \code{rr()} function (requires setting up an \eqn{n\times s} temporary data.frame). Scores and loadings form the matrices of these parameters within the CSAM.
#'
#' The returned list contains:
#' \itemize{
#'   \item \code{beta0}: species‑specific intercepts;
#'   \item \code{B}: a \eqn{g \times p} matrix of archetype slopes;
#'   \item \code{pi}: mixing proportions for the \eqn{g} archetypes;
#'   \item \code{U}: a \eqn{n \times d} matrix of site-specific factor scores;
#'   \item \code{Lambda}: a \eqn{s \times d} matrix of species-specific factor loadings;
#'   \item \code{phi}: species‑specific dispersion parameters (initially 1).
#' }
#'
#' An attribute \code{"sp clust"} stores the species‑level cluster assignments
#' from the k‑means step.
#'
#' @return A list containing initial parameter values for a SAM, with an
#'   attribute storing species cluster assignments.
#'
#' @export
#'
#' @importFrom glmmTMB glmmTMB getME
init.fa.pars <- function(Y, X, g = 3, family = poisson(), d = 2) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise usual SAM parameters
  start.pars = csam::init.sam.pars(Y, X, g = g, family = family)

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013
  beta0.init = vector("numeric", s)
  res = matrix(rep(0, n * s), n, s)
  for (sp in 1:s) {
    tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
    tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
    res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
  }

  # put data in long format
  tmp.dat = data.frame(res = as.vector(res), sp = factor(rep(1:s, each = nrow(res))), site = factor(rep(1:n, s)))
  # perform factor analysis on the residuals
  tmp.fa = glmmTMB::glmmTMB(res ~  0 + rr(sp + 0|site, d = d), data = tmp.dat, family = gaussian)

  # get the scores (these include zero's for full rank model)
  rr_b = glmmTMB::getME(tmp.fa, "b")
  # get out the observation-specific coefs (aka factor scores) by selecting d elements out of each s blocks
  scores = list()
  for (i in 1:d) {
    scores[[i]] = rr_b[seq(i,length(rr_b), by = s)]
  }
  start.pars$U = sapply(scores, cbind)
  start.pars$Lambda = tmp.fa$obj$env$report(tmp.fa$fit$parfull)$fact_load[[which(tmp.fa$modelInfo$grpVar == "site" & sapply(tmp.fa$modelInfo$reStruc$condReStruc, function(x){x$blockCode}) == 9)]]

  # return start pars in canonical order
  start.pars[c("beta0", "B", "pi", "U", "Lambda", "phi")]
}
