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
#' The penalized log-likelihood is
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
#'     and site scores \eqn{\mathbf{U}}, respectively, using penalized IRLS.
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
#' @param tol Numeric; convergence tolerance on the penalized log-likelihood.
#' @param verbose Logical; if `TRUE`, prints the penalized log-likelihood at
#'   each iteration.
#' @param init Optional list of initial values with elements
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
#'   \item{penalized_loglik}{Final penalized log-likelihood value.}
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
                 psi1 = 1, psi2 = 1, max_iter = 100, tol = 1e-6,
                 verbose = TRUE, init = NULL,
                 maxit_step1 = 1, maxit_step2 = 1, maxit_step3 = 1,
                 trace = TRUE) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    eps <- 1e-6
    if (family$family == "binomial") m <- pmin(pmax(m, eps), 1 - eps)
    if (family$family == "poisson") m <- pmax(m, eps)
    if (family$family == "gaussian") m <- m
    family$linkfun(m)
  }
  beta0 <- if (!is.null(init$beta0)) { init$beta0 } else { apply(Y, 2, safe_mean) }
  B <- if (!is.null(init$B)) { init$B } else { matrix(rnorm(g * p, 0, 0.1), g, p) }
  pi <- if (!is.null(init$pi)) { init$pi } else { rep(1 / g, g) }
  U <- if (!is.null(init$U)) { init$U } else { matrix(rnorm(n * d, 0, 0.1), n, d) }
  Lambda <- if (!is.null(init$Lambda)) { init$Lambda } else { matrix(rnorm(s * d, 0, 0.1), s, d) }
  phi <- if (!is.null(init$phi)) { init$phi } else { rep(1, s) }

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

  par.list = list(
    beta0 = beta0,
    B = B,
    pi = pi,
    U = U,
    Lambda = Lambda,
    phi = phi
  )

  prev_pll <- mzll(Y, X, par.list, g = g, d = d, family = family, psi1 = psi1, psi2 = psi2)
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

    if (iter < 1) {
      tau <- estep_post_probs(Y, X, par.list = par.list, family,
                              trunc.interval = TRUE)
    } else {
      tau <- estep_post_probs(Y, X, par.list = par.list, family)
    }
    # tau <- estep_post_probs(Y = Y, X = X, par.list = par.list, family = family)

    par.list$B <- mstep_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                         family = family, maxit = maxit_step1)

    sp <- mstep_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                             family = family, psi2 = psi2, maxit = maxit_step2)
    par.list$beta0  <- sp$beta0
    par.list$Lambda <- sp$Lambda
    par.list$phi    <- sp$phi

    par.list$U <- mstep_site_scores(Y = Y, X = X, par.list = par.list, tau = tau,
                           family = family, psi1 = psi1, maxit = maxit_step3)

    par.list$pi <- colMeans(tau)
    par.list$pi <- par.list$pi / sum(par.list$pi)

    # ll <- 0
    # for (j in 1:s) {
    #   yj <- Y[, j]
    #   lambda_j <- Lambda[j, , drop = FALSE]
    #   offset_common <- as.vector(U %*% t(lambda_j))
    #
    #   comp <- numeric(g)
    #   for (k in 1:g) {
    #     eta <- beta0[j] + X %*% B[k, ] + offset_common
    #     mu  <- family$linkinv(eta)
    #     dev <- family$dev.resids(y = yj, mu = mu, wt = rep(1, n))
    #     ll_comp <- -0.5 * family$aic(yj, rep(1, n), mu, rep(1, n), dev)
    #     comp[k] <- sum(ll_comp) + log(pi[k])
    #   }
    #   m <- max(comp)
    #   ll <- ll + m + log(sum(exp(comp - m)))
    # }
    #
    # pll <- ll - 0.5 * psi1 * sum(U^2) - 0.5 * psi2 * sum(Lambda^2)
    pll <-  mzll(Y, X, par.list, g = g, d = d, family = family, psi1 = psi1, psi2 = psi2)

    if (!is.finite(pll)) {
      warning("penalized log-likelihood became non-finite; stopping")
      break
    }

    if (verbose) {
      cat(sprintf("Iter %3d: penalized log-lik = %.6f\n", iter, pll))
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

      if (iter > 1) {
        diff <- pll - prev_pll
        trace_store$pll_diff[iter] <- diff
        trace_store$is_monotone[iter] <- (diff >= -1e-12)
      } else {
        trace_store$is_monotone[iter] <- TRUE
      }
    }

    if (abs(pll - prev_pll) < tol) break
    prev_pll <- pll
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

  out <- list(
    beta0 = par.list$beta0,
    B = par.list$B,
    pi = par.list$pi,
    U = par.list$U,
    Lambda = par.list$Lambda,
    phi = par.list$phi,
    tau = tau,
    par.list = par.list,
    penalized_loglik = pll,
    iterations = iter,
    family = family,
    Y = Y,
    X = X,
    psi1 = psi1,
    psi2 = psi2,
    g = g,
    d = d
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
#' @param tol Numeric; convergence tolerance on the penalized log-likelihood.
#' @param verbose Logical; if `TRUE`, prints the penalized log-likelihood at
#'   each iteration.
#' @param init Optional list of initial values with elements
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
#'   \item{penalized_loglik}{Final log-likelihood value. NOTE this is labelled "penalized" only for software compatibility with downstream functions for `csam`}
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
                verbose = TRUE, init = NULL,
                maxit_step1 = 1, maxit_step2 = 1,
                first_maxit_step1 = 50, first_maxit_step2 = 50,
                trace = TRUE) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    eps <- 1e-6
    if (family$family == "binomial") m <- pmin(pmax(m, eps), 1 - eps)
    if (family$family == "poisson") m <- pmax(m, eps)
    if (family$family == "gaussian") m <- m
    family$linkfun(m)
  }
  beta0 <- if (!is.null(init$beta0)) { init$beta0 } else { apply(Y, 2, safe_mean) }
  B <- if (!is.null(init$B)) { init$B } else { matrix(rnorm(g * p, 0, 0.1), g, p) }
  pi <- if (!is.null(init$pi)) { init$pi } else { rep(1 / g, g) }
  phi <- if (!is.null(init$phi)) { init$phi } else { rep(1, s) }

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

    if (iter == 1) {
      par.list$B <- mstep0_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                                     family = family, maxit = first_maxit_step1)
    } else {
      par.list$B <- mstep0_arch_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                                     family = family, maxit = maxit_step1)
    }

    if (iter == 1) {
      sp <- mstep0_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                                family = family, maxit = first_maxit_step2)
    } else {
      sp <- mstep0_species_pars(Y = Y, X = X, par.list = par.list, tau = tau,
                                family = family, maxit = maxit_step2)
    }

    par.list$beta0  <- sp$beta0
    par.list$phi    <- sp$phi

    par.list$pi <- colMeans(tau)
    par.list$pi <- par.list$pi / sum(par.list$pi)

    # ll <- 0
    # for (j in 1:s) {
    #   yj <- Y[, j]
    #
    #   comp <- numeric(g)
    #   for (k in 1:g) {
    #     eta <- beta0[j] + X %*% B[k, ]
    #     mu  <- family$linkinv(eta)
    #     dev <- family$dev.resids(y = yj, mu = mu, wt = rep(1, n))
    #     ll_comp <- -0.5 * family$aic(yj, rep(1, n), mu, rep(1, n), dev)
    #     comp[k] <- sum(ll_comp) + log(pi[k])
    #   }
    #   m <- max(comp)
    #   ll <- ll + m + log(sum(exp(comp - m)))
    # }
    #
    # pll <- ll
    pll <- mzll(Y, X, par.list, g = g, d = 1, family = family)

    if (verbose) {
      cat(sprintf("Iter %3d: log-lik = %.6f\n", iter, pll))
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

    if (abs(pll - prev_pll) < tol) { break }

    prev_pll <- pll

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
    penalized_loglik = pll,
    iterations = iter,
    family = family,
    Y = Y,
    X = X,
    g = g
  )

  if (trace) out$trace <- trace_store

  class(out) <- "csam"
  out
}
