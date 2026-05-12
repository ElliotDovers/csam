#' Fit a Correlated Species Archetype Model (CSAM)
#'
#' @description
#' Fits the correlated species archetype model
#' \deqn{
#'   h(\mu_{ijk}) = \beta_{0j} + \mathbf{x}_i^\top \boldsymbol{\beta}_k
#'                  + \mathbf{u}_i^\top \boldsymbol{\lambda}_j
#' }
#' for species \eqn{j = 1,\ldots,s}, sites \eqn{i = 1,\ldots,n}, and archetypes
#' \eqn{k = 1,\ldots,g}, where \eqn{h} is the GLM link function. The data are
#' modeled as a finite mixture of regression models with component densities
#' from a common exponential family.
#'
#' The penalised log-likelihood is
#' \deqn{
#'   \ell_{\text{pen}}(\theta; Y, X)
#'   = \ell(\theta; Y, X)
#'     - \frac{\psi_1}{2} \|\mathbf{U}\|_F^2
#'     - \frac{\psi_2}{2} \|\boldsymbol{\Lambda}\|_F^2,
#' }
#' where \eqn{\mathbf{U}} are site scores and \eqn{\boldsymbol{\Lambda}} are
#' species loadings, and \eqn{\|\cdot\|_F} is the Frobenius norm.
#'
#' Estimation proceeds via an ECM algorithm with:
#' \itemize{
#'   \item an E-step updating posterior species–archetype memberships
#'     \eqn{\tau_{jk}},
#'   \item three conditional M-steps updating archetype coefficients
#'     \eqn{\boldsymbol{B}}, species-specific parameters
#'     \eqn{(\boldsymbol{\beta}_0, \boldsymbol{\Lambda}, \boldsymbol{\phi})},
#'     and site scores \eqn{\mathbf{U}}, respectively, using penalised IRLS.
#' }
#'
#' @param Y A numeric matrix of responses of dimension \eqn{n \times s},
#'   with rows corresponding to sites and columns to species.
#' @param X A numeric matrix of predictors of dimension \eqn{n \times p}.
#' @param g Integer; number of archetypes (mixture components).
#' @param d Integer; dimension of the common factor-analytic term.
#' @param family A GLM family object (e.g. [stats::poisson()], [stats::binomial()],
#'   [stats::gaussian()]) specifying the exponential family and link.
#' @param psi1 Non-negative numeric; ridge penalty parameter on the site scores
#'   matrix \eqn{\mathbf{U}}.
#' @param psi2 Non-negative numeric; ridge penalty parameter on the loadings
#'   matrix \eqn{\boldsymbol{\Lambda}}.
#' @param max_iter Integer; maximum number of ECM iterations.
#' @param tol Numeric; convergence tolerance on the penalised log-likelihood.
#' @param verbose Logical; if `TRUE`, prints the penalised log-likelihood at
#'   each iteration.
#' @param start Optional list of initial values with elements
#'   `beta0`, `B`, `pi`, `U`, `Lambda`, `phi`. Any missing elements are
#'   initialized internally.
#' @param maxit_step1 Integer; maximum P-IRLS iterations in the archetype
#'   update step.
#' @param maxit_step2 Integer; maximum P-IRLS iterations in the species
#'   parameter update step.
#' @param maxit_step3 Integer; maximum P-IRLS iterations in the site scores
#'   update step.
#' @param trace Logical; if `TRUE`, stores parameter and likelihood traces
#'   across ECM iterations for diagnostics and plotting.
#' @param backend Character string; if `TRUE`, stores parameter and likelihood traces
#'   across ECM iterations for diagnostics and plotting.
#' @param constrain Logical; if `TRUE`, performs an identifiability constraint on the factor analytic term.
#' @param inner.constrain Logical; if `TRUE`, performs an identifiability constraint on the factor analytic term within each iteration of the EM algorithm.
#'
#' @return An object of class `"csam"` with components:
#' \describe{
#'   \item{beta0}{Length-\(s\) vector of species intercepts.}
#'   \item{B}{\eqn{g \times p} matrix of archetype-specific regression
#'     coefficients.}
#'   \item{pi}{Length-\(g\) vector of mixing proportions.}
#'   \item{U}{\eqn{n \times d} matrix of site scores.}
#'   \item{Lambda}{\eqn{s \times d} matrix of species loadings.}
#'   \item{phi}{Length-\(s\) vector of species-specific dispersion parameters
#'     (for Gaussian family).}
#'   \item{tau}{\eqn{s \times g} matrix of posterior species–archetype
#'     memberships.}
#'   \item{penalised_loglik}{Final penalised log-likelihood value.}
#'   \item{iterations}{Number of ECM iterations performed.}
#'   \item{family}{The GLM family used.}
#'   \item{trace}{If `trace = TRUE`, a list containing parameter and
#'     likelihood traces across iterations.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 50; s <- 10; p <- 2
#' X <- cbind(1, scale(matrix(rnorm(n * (p - 1)), n, p - 1)))
#' Y <- matrix(rpois(n * s, lambda = 5), n, s)
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' plot.csam(fit, param = "pll")
#' }
#'
#' @export
csam <- function(Y, X, g = 3, d = 2, family = poisson(),
                 psi1 = 0, psi2 = 0, max_iter = 100, tol = 1e-3,
                 verbose = TRUE, start = NULL,
                 maxit_step1 = 5, maxit_step2 = 5, maxit_step3 = 5,
                 trace = TRUE, backend = c("C++", "R"), constrain = FALSE, inner.constrain = FALSE, starts.at.steps = FALSE, trunc.tau.until.iter = 2, project.loadings = FALSE) {

  backend <- match.arg(backend)
  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # check if the supplied starting parameters are a fitted csam object and correct as needed
  if (!is.null(start)) {
    if (class(start) == "csam") {
      start = start$par.list
    }
  }

  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    eps <- 1e-6
    if (family$family == "binomial") m <- pmin(pmax(m, eps), 1 - eps)
    if (family$family == "poisson") m <- pmax(m, eps)
    if (family$family == "gaussian") m <- m
    family$linkfun(m)
  }
  beta0 <- if (!is.null(start$beta0)) { start$beta0 } else { apply(Y, 2, safe_mean) }
  B <- if (!is.null(start$B)) { start$B } else { matrix(rnorm(g * p, 0, 0.1), g, p) }
  pi <- if (!is.null(start$pi)) { start$pi } else { rep(1 / g, g) }
  U <- if (!is.null(start$U)) { start$U } else { matrix(rnorm(n * d, 0, 0.1), n, d) }
  Lambda <- if (!is.null(start$Lambda)) { start$Lambda } else { matrix(rnorm(s * d, 0, 0.1), s, d) }
  phi <- if (!is.null(start$phi)) { start$phi } else { rep(1, s) }

  # some checks on warm starts
  if (g != nrow(B)) {
    stop("mismatch between initial archetype parameters, B, and specified g")
  }
  if (p != ncol(B)) {
    stop("mismatch between initial predictor parameters, B, and columns in supplied X")
  }
  if (s != length(beta0)) {
    stop("mismatch between initial response intercepts, beta0, and response columns in supplied Y")
  }
  if (s != length(phi)) {
    stop("mismatch between initial response dispersion parameters, phi, and response columns in supplied Y")
  }
  if (g != length(pi)) {
    stop("mismatch between initial mixing parameters, pi, and specified g")
  }
  if (s != nrow(Lambda)) {
    stop("mismatch between initial factor loadings, Lambda, and response columns in supplied Y")
  }
  if (n != nrow(U)) {
    stop("mismatch between initial factor scores, U, and response rows in supplied Y")
  }
  if (d != ncol(Lambda)) {
    stop("mismatch between initial factor loadings, Lambda, and specified d")
  }
  if (d != ncol(U)) {
    stop("mismatch between initial factor scores, U, and specified d")
  }

  prev_par.list <- par.list <- list(
    beta0 = beta0,
    B = B,
    pi = pi,
    U = U,
    Lambda = Lambda,
    phi = phi
  )

  prev_pll <- mzll(Y, X, par.list,
                   # list(beta0 = beta0, B = B, pi = pi, U = U, Lambda = Lambda, phi = phi),
                   g = g, d = d, family = family, psi1 = psi1, psi2 = psi2)
  if (verbose) {
    cat(sprintf("Iter %3d: log-lik = %.6f\n", 0, prev_pll))
  }

  if (trace) {
    trace_store <- list(
      beta0    = vector("list", max_iter),
      B        = vector("list", max_iter),
      pi       = vector("list", max_iter),
      U        = vector("list", max_iter),
      Lambda   = vector("list", max_iter),
      phi      = vector("list", max_iter),
      tau      = vector("list", max_iter),
      pll      = rep(NA_real_, max_iter),
      pll_diff = rep(NA_real_, max_iter),
      is_monotone = rep(NA, max_iter)
    )
  }

  for (iter in 1:max_iter) {

    if (iter < trunc.tau.until.iter) {
      tau <- estep_post_probs(Y, X, par.list = par.list, family,
                              trunc.interval = TRUE)
    } else {
      tau <- estep_post_probs(Y, X, par.list = par.list, family)
    }
    # tau <- estep_post_probs(Y = Y, X = X, par.list = par.list, family = family)

    par.list$B <- mstep_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                         family = family, maxit = maxit_step1, use.starts = starts.at.steps)

    ############################################################################
    # can we constrain Lambda to be orthogonal to the archetype parameters?
    if (project.loadings) {
      # After updating U and Lambda
      z_hat <- apply(tau, 1, which.max)

      # Species-level archetype effects
      B_sp <- par.list$B[z_hat, , drop = FALSE]
      b <- B_sp / sqrt(sum(B_sp^2))

      # Orthogonalise Lambda wrt archetype space
      Lambda <- Lambda - b %*% crossprod(b, Lambda)

      par.list$Lambda <- Lambda
    }
    ############################################################################

    sp <- mstep_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                             family = family, psi2 = psi2, maxit = maxit_step2, backend = backend, use.starts = starts.at.steps)
    par.list$beta0  <- sp$beta0
    par.list$Lambda <- sp$Lambda
    par.list$phi    <- sp$phi

    par.list$U <- mstep_site_scores(Y = Y, X = X, par.list = par.list, tau = tau,
                           family = family, psi1 = psi1, maxit = maxit_step3, backend = backend, use.starts = starts.at.steps)

    # apply constraints to the factor analytic terms if desired
    if (inner.constrain) {
      if (all(par.list$U == 0)) {
        warning("no constraint applied to factor terms: all scores are zero")
      } else {
        try(assign("new.fa", correct.uv(par.list$U, par.list$Lambda)))
        if (exists("new.fa")) {
          par.list$U <- new.fa$u
          par.list$Lambda <- new.fa$v
        } else {
          warning("inner constraint could not be applied")
        }
      }
    }

    par.list$pi <- colMeans(tau)
    par.list$pi <- par.list$pi / sum(par.list$pi)

    pll <-  mzll(Y, X, par.list,
                 # list(beta0 = beta0, B = B, pi = pi, U = U, Lambda = Lambda, phi = phi),
                 g = g, d = d, family = family, psi1 = psi1, psi2 = psi2)

    if (verbose) {
      cat(sprintf(paste0("Iter %3d: ", if (psi1 != 0 | psi2 != 0) {"penalised "} else {""}, "logLik = %.6f | Ratio = %3f\n"), iter, pll, exp(prev_pll - pll)))
    }

    if (trace) {
      trace_store$beta0[[iter]]  <- par.list$beta0
      trace_store$B[[iter]]      <- par.list$B
      trace_store$pi[[iter]]     <- par.list$pi
      trace_store$U[[iter]]      <- par.list$U
      trace_store$Lambda[[iter]] <- par.list$Lambda
      trace_store$phi[[iter]]    <- par.list$phi
      trace_store$tau[[iter]]    <- tau
      trace_store$pll[iter]      <- pll
      diff <- pll - prev_pll
      trace_store$pll_diff[iter] <- diff
      trace_store$is_monotone[iter] <- (diff >= -1e-12)
    }

    if (!is.finite(pll)) {
      warning(paste0(if (psi1 != 0 | psi2 != 0) {"penalised "} else {""}, "log-likelihood became non-finite; stopping at previous iteration"))
      conv = NA
      par.list <- prev_par.list
      pll <- prev_pll
      iter <- iter - 1
      break
    }

    # if (pll - prev_pll < -1e-12) {
    #   warning(paste0(if (psi1 != 0 | psi2 != 0) {"penalised "} else {""}, "log-likelihood no longer monotonic increasing; stopping at previous iteration"))
    #   conv = 2
    #   par.list <- prev_par.list
    #   pll <- prev_pll
    #   iter <- iter - 1
    #   break
    # }

    if (abs(pll - prev_pll) < tol) {
      conv = 0
      break
    }

    prev_pll <- pll
    prev_par.list <- par.list
    conv = 1
  }

  if (trace) {
    trace_store$beta0    <- trace_store$beta0[1:iter]
    trace_store$B        <- trace_store$B[1:iter]
    trace_store$pi       <- trace_store$pi[1:iter]
    trace_store$U        <- trace_store$U[1:iter]
    trace_store$Lambda   <- trace_store$Lambda[1:iter]
    trace_store$phi      <- trace_store$phi[1:iter]
    trace_store$tau      <- trace_store$tau[1:iter]
    trace_store$pll      <- trace_store$pll[1:iter]
    trace_store$pll_diff <- trace_store$pll_diff[1:iter]
    trace_store$is_monotone <- trace_store$is_monotone[1:iter]
  }

  # apply constraints to the factor analytic terms if desired
  if (constrain) {
    if (all(par.list$U == 0)) {
      warning("no constraint applied to factor terms: all scores are zero")
    } else {
      try(assign("new.fa", correct.uv(par.list$U, par.list$Lambda)))
      if (exists("new.fa")) {
        par.list$U <- new.fa$u
        par.list$Lambda <- new.fa$v
      } else {
        warning("post-hoc constraint could not be applied")
      }
    }
  }

  out <- list(
    beta0 = par.list$beta0,
    B = par.list$B,
    pi = par.list$pi,
    U = par.list$U,
    Lambda = par.list$Lambda,
    phi = par.list$phi,
    tau = tau,
    par.list = par.list,
    penalised_loglik = pll,
    marginalised_loglik = NA,
    iterations = iter,
    family = family,
    Y = Y,
    X = X,
    psi1 = psi1,
    psi2 = psi2,
    g = g,
    d = d,
    convergence = conv
  )

  if (trace) out$trace <- trace_store

  class(out) <- "csam"
  out
}

#' Fit a vanilla Species Archetype Model (SAM)
#'
#' @description
#' Fits the species archetype model
#' \deqn{
#'   h(\mu_{ijk}) = \beta_{0j} + \mathbf{x}_i^\top \boldsymbol{\beta}_k
#' }
#' for species \eqn{j = 1,\ldots,s}, sites \eqn{i = 1,\ldots,n}, and archetypes
#' \eqn{k = 1,\ldots,g}, where \eqn{h} is the GLM link function. The data are
#' modeled as a finite mixture of regression models with component densities
#' from a common exponential family.
#'
#' The log-likelihood is
#' \deqn{
#'   \ell\left(\boldsymbol{\theta};\boldsymbol{Y},\boldsymbol{X}\right) =\sum_{j=1}^{s}\log\left\{\sum_{k=1}^{g}\pi_{k}\prod_{i=1}^{n}f\left(y_{ij}|\mu_{ijk},\phi_{j}\right)\right\}
#' }
#'
#' Estimation proceeds via an ECM algorithm with:
#' \itemize{
#'   \item an E-step updating posterior species–archetype memberships
#'     \eqn{\tau_{jk}},
#'   \item two conditional M-steps updating archetype coefficients
#'     \eqn{\boldsymbol{B}}, and species-specific parameters
#'     \eqn{(\boldsymbol{\beta}_0, \boldsymbol{\phi})} using IRLS.
#' }
#'
#' @param Y A numeric matrix of responses of dimension \eqn{n \times s},
#'   with rows corresponding to sites and columns to species.
#' @param X A numeric matrix of predictors of dimension \eqn{n \times p}.
#' @param g Integer; number of archetypes (mixture components).
#' @param family A GLM family object (e.g. [stats::poisson()], [stats::binomial()],
#'   [stats::gaussian()]) specifying the exponential family and link.
#' @param max_iter Integer; maximum number of ECM iterations.
#' @param tol Numeric; convergence tolerance on the penalised log-likelihood.
#' @param verbose Logical; if `TRUE`, prints the penalised log-likelihood at
#'   each iteration.
#' @param start Optional list of initial values with elements
#'   `beta0`, `B`, `pi`, `phi`. Any missing elements are
#'   initialized internally.
#' @param maxit_step1 Integer; maximum P-IRLS iterations in the archetype
#'   update step.
#' @param maxit_step2 Integer; maximum P-IRLS iterations in the species
#'   parameter update step.
#' @param trace Logical; if `TRUE`, stores parameter and likelihood traces
#'   across ECM iterations for diagnostics and plotting.
#'
#' @return An object of class `"csam"` with components:
#' \describe{
#'   \item{beta0}{Length-\(s\) vector of species intercepts.}
#'   \item{B}{\eqn{g \times p} matrix of archetype-specific regression
#'     coefficients.}
#'   \item{pi}{Length-\(g\) vector of mixing proportions.}
#'   \item{phi}{Length-\(s\) vector of species-specific dispersion parameters
#'     (for Gaussian family).}
#'   \item{tau}{\eqn{s \times g} matrix of posterior species–archetype
#'     memberships.}
#'   \item{penalised_loglik}{Final log-likelihood value. NOTE this is labelled "penalised" only for software compatibility with downstream functions for `csam`}
#'   \item{iterations}{Number of ECM iterations performed.}
#'   \item{family}{The GLM family used.}
#'   \item{trace}{If `trace = TRUE`, a list containing parameter and
#'     likelihood traces across iterations.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 50; s <- 10; p <- 2
#' X <- cbind(1, scale(matrix(rnorm(n * (p - 1)), n, p - 1)))
#' Y <- matrix(rpois(n * s, lambda = 5), n, s)
#' fit <- sam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' plot.csam(fit, param = "ll")
#' }
#'
#' @export
sam <- function(Y, X, g = 3, family = poisson(),
                max_iter = 100, tol = 1e-6,
                verbose = TRUE, start = NULL,
                maxit_step1 = 5, maxit_step2 = 5,
                #first_maxit_step1 = 50, first_maxit_step2 = 50,
                trace = TRUE) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # check if the supplied starting parameters are a fitted csam object and correct as needed
  if (!is.null(start)) {
    if (class(start) == "csam") {
      start = start$par.list
    }
  }

  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    eps <- 1e-6
    if (family$family == "binomial") m <- pmin(pmax(m, eps), 1 - eps)
    if (family$family == "poisson") m <- pmax(m, eps)
    if (family$family == "gaussian") m <- m
    family$linkfun(m)
  }
  beta0 <- if (!is.null(start$beta0)) { start$beta0 } else { apply(Y, 2, safe_mean) }
  B <- if (!is.null(start$B)) { start$B } else { matrix(rnorm(g * p, 0, 0.1), g, p) }
  pi <- if (!is.null(start$pi)) { start$pi } else { rep(1 / g, g) }
  phi <- if (!is.null(start$phi)) { start$phi } else { rep(1, s) }

  # some checks on warm starts
  if (g != nrow(B)) {
    stop("mismatch between initial archetype parameters, B, and specified g")
  }
  if (p != ncol(B)) {
    stop("mismatch between initial predictor parameters, B, and columns in supplied X")
  }
  if (s != length(beta0)) {
    stop("mismatch between initial response intercepts, beta0, and response columns in supplied Y")
  }
  if (s != length(phi)) {
    stop("mismatch between initial response dispersion parameters, phi, and response columns in supplied Y")
  }
  if (g != length(pi)) {
    stop("mismatch between initial mixing parameters, pi, and specified g")
  }

  par.list = list(
    beta0 = beta0,
    B = B,
    pi = pi,
    phi = phi
  )

  prev_pll <- mzll(Y, X, par.list, g = g, d = 1, family = family)
  if (verbose) {
    cat(sprintf("Iter %3d: log-lik = %.6f\n", 0, prev_pll))
  }


  if (trace) {
    trace_store <- list(
      beta0    = vector("list", max_iter),
      B        = vector("list", max_iter),
      pi       = vector("list", max_iter),
      phi      = vector("list", max_iter),
      tau      = vector("list", max_iter),
      pll      = rep(NA_real_, max_iter),
      pll_diff = rep(NA_real_, max_iter),
      is_monotone = rep(NA, max_iter)
    )
  }

  for (iter in 1:max_iter) {

    if (iter == 1) {
      tau <- estep0_post_probs(Y, X, par.list = par.list, family,
                              trunc.interval = TRUE)
    } else {
      tau <- estep0_post_probs(Y, X, par.list = par.list, family)
    }
    # tau <- estep0_post_probs(Y = Y, X = X, par.list = par.list, family = family)

    # if (iter == 1) {
    #   par.list$B <- mstep0_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
    #                                  family = family, maxit = first_maxit_step1)
    # } else {
    #   par.list$B <- mstep0_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
    #                                  family = family, maxit = maxit_step1)
    # }
    par.list$B <- mstep0_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                                   family = family, maxit = maxit_step1)

    # if (iter == 1) {
    #   sp <- mstep0_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
    #                             family = family, maxit = first_maxit_step2)
    # } else {
    #   sp <- mstep0_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
    #                             family = family, maxit = maxit_step2)
    # }
    sp <- mstep0_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                              family = family, maxit = maxit_step2)

    par.list$beta0  <- sp$beta0
    par.list$phi    <- sp$phi

    par.list$pi <- colMeans(tau)
    par.list$pi <- par.list$pi / sum(par.list$pi)

    pll <- mzll(Y, X, par.list, g = g, d = 1, family = family)

    if (!is.finite(pll)) {
      warning("log-likelihood became non-finite; stopping")
      conv = NA
      break
    }

    if (verbose) {
      cat(sprintf("Iter %3d: log-lik = %.6f | Ratio = %3f\n", iter, pll, exp(prev_pll - pll)))
    }

    if (trace) {
      trace_store$beta0[[iter]]  <- par.list$beta0
      trace_store$B[[iter]]      <- par.list$B
      trace_store$pi[[iter]]     <- par.list$pi
      trace_store$phi[[iter]]    <- par.list$phi
      trace_store$tau[[iter]]    <- tau
      trace_store$pll[iter]      <- pll

      if (iter > 1) {
        diff <- pll - prev_pll
        trace_store$pll_diff[iter] <- diff
        trace_store$is_monotone[iter] <- (diff >= -1e-12)
      } else {
        trace_store$is_monotone[iter] <- TRUE
      }
    }

    if (abs(pll - prev_pll) < tol) {
      conv = 0
      break
    }

    prev_pll <- pll
    conv = 1

  }

  if (trace) {
    trace_store$beta0    <- trace_store$beta0[1:iter]
    trace_store$B        <- trace_store$B[1:iter]
    trace_store$pi       <- trace_store$pi[1:iter]
    trace_store$phi      <- trace_store$phi[1:iter]
    trace_store$tau      <- trace_store$tau[1:iter]
    trace_store$pll      <- trace_store$pll[1:iter]
    trace_store$pll_diff <- trace_store$pll_diff[1:iter]
    trace_store$is_monotone <- trace_store$is_monotone[1:iter]
  }

  out <- list(
    beta0 = par.list$beta0,
    B = par.list$B,
    pi = par.list$pi,
    phi = par.list$phi,
    tau = tau,
    par.list = par.list,
    penalised_loglik = pll,
    iterations = iter,
    family = family,
    Y = Y,
    X = X,
    g = g,
    convergence = conv
  )

  if (trace) out$trace <- trace_store

  class(out) <- "csam"
  out
}

#' Fit a Correlated Species Archetype Model (CSAM)
#'
#' @description
#' Fits the correlated species archetype model
#' \deqn{
#'   h(\mu_{ijk}) = \beta_{0j} + \mathbf{x}_i^\top \boldsymbol{\beta}_k
#'                  + \mathbf{u}_i^\top \boldsymbol{\lambda}_j
#' }
#' for species \eqn{j = 1,\ldots,s}, sites \eqn{i = 1,\ldots,n}, and archetypes
#' \eqn{k = 1,\ldots,g}, where \eqn{h} is the GLM link function. The data are
#' modeled as a finite mixture of regression models with component densities
#' from a common exponential family.
#'
#' Estimation proceeds via a gradient-based optiimiser with:
#'
#' The penalised log-likelihood is
#' \deqn{
#'   \ell_{\text{pen}}(\theta; Y, X)
#'   = \ell(\theta; Y, X)
#'     - \frac{\psi_1}{2} \|\mathbf{U}\|_F^2
#'     - \frac{\psi_2}{2} \|\boldsymbol{\Lambda}\|_F^2,
#' }
#' where \eqn{\mathbf{U}} are site scores and \eqn{\boldsymbol{\Lambda}} are
#' species loadings, and \eqn{\|\cdot\|_F} is the Frobenius norm. Alternatively,
#' the \eqn{\boldsymbol{U}} may be considered as standard normal variables and the
#' model fitted via a marginalised log-likelihood approximated by Laplace (as in TMB)
#'
#' @param Y A numeric matrix of responses of dimension \eqn{n \times s},
#'   with rows corresponding to sites and columns to species.
#' @param X A numeric matrix of predictors of dimension \eqn{n \times p}.
#' @param g Integer; number of archetypes (mixture components).
#' @param d Integer; dimension of the common factor-analytic term.
#' @param family A GLM family object (e.g. [stats::poisson()], [stats::binomial()],
#'   [stats::gaussian()]) specifying the exponential family and link.
#' @param psi1 Non-negative numeric; ridge penalty parameter on the site scores
#'   matrix \eqn{\mathbf{U}}.
#' @param psi2 Non-negative numeric; ridge penalty parameter on the loadings
#'   matrix \eqn{\boldsymbol{\Lambda}}.
#' @param U.random Logical; indicating whether to treat the site effects, \eqn{\boldsymbol{U}} as random or fixed effects.
#' @param se Logical; if `TRUE`, approximates standard errors using [TMB::sdreport()]
#' @param nlminb.control list; of control parameters for [stats::nlminb()]
#' @param start Optional list of initial values with elements
#'   `beta0`, `B`, `pi`, `U`, `Lambda`, `phi`. Any missing elements are
#'   initialized internally. Alternatively, a fitted \code{csam} model object
#' @param vanilla.sam Logical; if `TRUE`, the a standard SAM will be fitted that excludes residual correlation (i.e. \eqn{\boldsymbol{U}\boldsymbol{\Lambda}^{\top}})
#'
#' @return An object of class `"csam"` with components:
#' \describe{
#'   \item{beta0}{Length-\(s\) vector of species intercepts.}
#'   \item{B}{\eqn{g \times p} matrix of archetype-specific regression
#'     coefficients.}
#'   \item{pi}{Length-\(g\) vector of mixing proportions.}
#'   \item{U}{\eqn{n \times d} matrix of site scores.}
#'   \item{Lambda}{\eqn{s \times d} matrix of species loadings.}
#'   \item{phi}{Length-\(s\) vector of species-specific dispersion parameters
#'     (for Gaussian family).}
#'   \item{tau}{not available}
#'   \item{penalised_loglik}{Final penalised log-likelihood value.}
#'   \item{iterations}{Number of ECM iterations performed.}
#'   \item{family}{The GLM family used.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 50; s <- 10; p <- 2
#' X <- cbind(1, scale(matrix(rnorm(n * (p - 1)), n, p - 1)))
#' Y <- matrix(rpois(n * s, lambda = 5), n, s)
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' plot.csam(fit, param = "pll")
#' }
#'
#' @export
#'
#' @importFrom stats binomial nlminb
#' @importFrom TMB MakeADFun sdreport
csam_tmb <- function(Y, X, g = 3, d = 2, family = poisson(),
                 psi1 = 1, psi2 = 1, U.random = FALSE, se = FALSE,
                 nlminb.control = list(eval.max = 200, iter.max = 150, trace = FALSE),
                 start = NULL, vanilla.sam = FALSE) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # check if the supplied starting parameters are a fitted csam object and correct as needed
  if (!is.null(start)) {
    if (class(start) == "csam") {
      start = start$par.list
    }
  }

  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    eps <- 1e-6
    if (family$family == "binomial") m <- pmin(pmax(m, eps), 1 - eps)
    if (family$family == "poisson") m <- pmax(m, eps)
    if (family$family == "gaussian") m <- m
    family$linkfun(m)
  }
  beta0 <- if (!is.null(start$beta0)) { start$beta0 } else { apply(Y, 2, safe_mean) }
  B <- if (!is.null(start$B)) { start$B } else { matrix(rnorm(g * p, 0, 0.1), g, p) }
  pi <- if (!is.null(start$pi)) { start$pi } else { rep(1 / g, g) }
  U <- if (!is.null(start$U)) { start$U } else { matrix(rnorm(n * d, 0, 0.1), n, d) }
  Lambda <- if (!is.null(start$Lambda)) { start$Lambda } else { matrix(rnorm(s * d, 0, 0.1), s, d) }
  phi <- if (!is.null(start$phi)) { start$phi } else { rep(1, s) }

  # adjust the start parameters for factor analytic terms when fitting a vanilla SAM
  if (vanilla.sam) {
    # set both to zero which will be mapped in the objective function later
    Lambda = matrix(rep(0, s * d), s, d)
    U <- matrix(rep(0, n * d), n, d)
  }

  # some checks on warm starts
  if (g != nrow(B)) {
    stop("mismatch between initial archetype parameters, B, and specified g")
  }
  if (p != ncol(B)) {
    stop("mismatch between initial predictor parameters, B, and columns in supplied X")
  }
  if (s != length(beta0)) {
    stop("mismatch between initial response intercepts, beta0, and response columns in supplied Y")
  }
  if (s != length(phi)) {
    stop("mismatch between initial response dispersion parameters, phi, and response columns in supplied Y")
  }
  if (g != length(pi)) {
    stop("mismatch between initial mixing parameters, pi, and specified g")
  }
  if (s != nrow(Lambda)) {
    stop("mismatch between initial factor loadings, Lambda, and response columns in supplied Y")
  }
  if (n != nrow(U)) {
    stop("mismatch between initial factor scores, U, and response rows in supplied Y")
  }
  if (d != ncol(Lambda)) {
    stop("mismatch between initial factor loadings, Lambda, and specified d")
  }
  if (d != ncol(U)) {
    stop("mismatch between initial factor scores, U, and specified d")
  }

  # set up the appropriate TMB lists for the objective function
  par.list.init = list(
    beta0 = beta0,
    B = B,
    logit_pi = stats::binomial()$linkfun(pi[1:(g - 1)]),
    U = U,
    Lambda = Lambda,
    log_phi = log(phi)
  )
  data.list <- list(
    Y = Y,
    X =  X,
    psi1 = if (!vanilla.sam) { psi1 } else { 0 },
    psi2 = if (!vanilla.sam) { psi2 } else { 0 },
    family = switch(family$family,
                    poisson = 0,
                    binomial = 1,
                    gaussian = 2),  # 0=Poisson, 1=Binomial, 2=Gaussian,...
    lik_type = as.numeric(U.random) # 0 = penalised, 1 = marginalised
  )

  # create the mapping of parameters according to csam type and family
  mapping = list()
  if (family$family %in% c("binomial", "poisson")) {
    mapping$log_phi = factor(rep(NA, s))
  }
  if (vanilla.sam) {
    mapping$U = factor(rep(NA, n * d))
    mapping$Lambda = factor(rep(NA, s * d))
  }

  # set up the objective function
  if (U.random) {
    if (vanilla.sam) {
      warning("Marginalised likelihood for a vanilla SAM will be slightly off")
    }
    obj <- TMB::MakeADFun(
      data = data.list,
      parameters = par.list.init,
      DLL = "csam", random = "U", map = mapping,
      silent = TRUE
    )
  } else {
    obj <- TMB::MakeADFun(
      data = data.list,
      parameters = par.list.init,
      DLL = "csam", map = mapping,
      silent = TRUE
    )
  }

  # fit the model
  opt = stats::nlminb(obj$par, obj$fn, obj$gr, control = nlminb.control)

  # adjust the parameters to the canonical setting
  par.list = lapply(split(opt$par, names(opt$par)), unname)
  par.list$B = matrix(par.list$B, nrow = g)
  if (U.random) {
    par.list$U = obj$env$parList()$U
  }
  if (!is.null(par.list$U)) {
    par.list$U = matrix(par.list$U, nrow = n)
  }
  if (!is.null(par.list$Lambda)) {
    par.list$Lambda = matrix(par.list$Lambda, nrow = s)
  }
  par.list$pi = binomial()$linkinv(par.list$logit_pi)
  par.list$pi[g] = 1 - sum(par.list$pi[1:(g - 1)])
  par.list$logit_pi = NULL
  if (family$family %in% c("binomial", "poisson")) {
    par.list$phi = rep(1, s)
  } else {
    par.list$phi = exp(par.list$log_phi)
    par.list$log_phi = NULL
  }

  # set up the output list to match other csam model objects
  out = par.list
  out$tau = tau_from_tmb_fit(Y, X, par.list, family)
  out$par.list = par.list
  out$penalised_loglik = if (U.random) {NA} else {-opt$objective}
  out$marginalised_loglik = if (U.random) {-opt$objective} else {NA}
  out$iterations = opt$iterations
  out$family = family
  out$Y = Y
  out$X = X
  out$psi1 = psi1
  out$psi2 = psi2
  out$g = g
  out$d = d
  out$convergence = opt$convergence

  # get standard errors if required
  if (se) {
    rep = TMB::sdreport(obj)
    out$se = summary(rep)
  }

  class(out) <- "csam"
  out
}
