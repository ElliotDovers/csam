#' Fit a penalised generalised linear model by penalised IRLS
#'
#' Fit a generalised linear model with a quadratic penalty on coefficients
#' using an iteratively reweighted least squares (IRLS) approach. The penalty
#' is specified via a diagonal penalty weight vector and a scalar tuning
#' parameter \code{lambda}; the effective penalty matrix is
#' \eqn{\sqrt{\lambda}\,P} where \code{P} = diag(\code{penalty_weights}).
#'
#' @param X **matrix** of predictors with \eqn{n} rows (observations) and
#'   \eqn{p} columns (features).
#' @param y **numeric** response vector of length \eqn{n}.
#' @param weights **numeric** optional observation weights of length \eqn{n}.
#'   Defaults to \code{rep(1, n)}.
#' @param offset **numeric** optional offset vector of length \eqn{n}.
#'   Defaults to \code{rep(0, n)}.
#' @param penalty_weights **numeric** nonnegative vector of length \eqn{p}
#'   giving diagonal entries of the penalty matrix \code{P}. Defaults to
#'   \code{rep(0, p)} (no penalty).
#' @param lambda **numeric** nonnegative scalar that scales the penalty matrix.
#'   Defaults to \code{0} (no penalty).
#' @param family a \pkg{stats}-style family object (for example
#'   \code{gaussian()}, \code{binomial()}, \code{poisson()}). The function
#'   uses \code{family$linkinv}, \code{family$mu.eta} and
#'   \code{family$variance}. Default is \code{gaussian()}.
#' @param start **numeric** optional starting values for the coefficients of
#'   length \eqn{p}. If \code{NULL} (default) coefficients are initialized to
#'   zero. A length mismatch raises an error.
#' @param maxit **integer** maximum number of IRLS iterations. Default is
#'   \code{50}.
#' @param tol **numeric** convergence tolerance on the maximum absolute change
#'   in coefficients between iterations. Default is \code{1e-6}.
#' @param backend **character** one of \code{"C++"} or \code{"R"} indicating the framework in which the P-IRLS scheme is executed.
#'   Default is \code{"C++"} which is more efficient, \code{"R"} is only retained for legacy and checking.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{coefficients}}{numeric vector of length \eqn{p} with fitted
#'     regression coefficients.}
#'   \item{\code{converged}}{logical indicating whether the algorithm converged.}
#'   \item{\code{iterations}}{numeric indicating the number of iterations used.}
#' }
#'
#' @details
#' The algorithm performs P-IRLS where, at each iteration, a weighted least
#' squares problem with augmented rows is solved:
#' \preformatted{
#'   A = rbind(sqrt(W) * X, sqrt(lambda) * P)
#'   b = c(sqrt(W) * z, rep(0, p))
#' }
#' where \code{W} are the IRLS weights and \code{z} is the working response.
#' The penalty matrix \code{P} is diagonal with entries given by
#' \code{penalty_weights}. When \code{lambda == 0} or all
#' \code{penalty_weights == 0} the fit reduces to an unpenalised IRLS fit.
#'
#' If \code{start} is provided it must have length equal to the number of
#' columns of \code{X}; otherwise an error is raised.
#'
#' @examples
#' ## Gaussian example
#' set.seed(1)
#' n <- 100; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(1, -0.5, 0, 0, 0)
#' y <- X %*% beta_true + rnorm(n)
#' fit <- penalised_glm_fit(X, y, lambda = 1, penalty_weights = c(0, 1, 1, 1, 1))
#' print(fit$coefficients)
#'
#' ## Binomial example (logit link)
#' set.seed(2)
#' n <- 200; p <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(0.5, -1, 0.2)
#' eta <- X %*% beta_true
#' p_true <- 1 / (1 + exp(-eta))
#' y_bin <- rbinom(n, size = 1, prob = p_true)
#' fit_bin <- penalised_glm_fit_R(X, y_bin, family = binomial(), lambda = 0.5,
#'                              penalty_weights = rep(1, p))
#' head(fit_bin$mu)
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#' @export
penalised_glm_fit <- function(X, y, weights = NULL, offset = NULL,
                              penalty_weights = NULL, lambda = 0,
                              family = gaussian(), start = NULL,
                              maxit = 50, tol = 1e-6, verbose = FALSE, backend = c("C++", "R")) {
  backend <- match.arg(backend)
  if (backend == "R") {
    out <- penalised_glm_fit_R(X = X, y = y, weights = weights, offset = offset,
                        penalty_weights = penalty_weights, lambda = lambda,
                        family = family, start = start,
                        maxit = maxit, tol = tol, verbose = verbose)
  }
  else {
    out <- penalised_glm_fit_C(X = X, y = y, weights = weights, offset = offset,
                               penalty_weights = penalty_weights, lambda = lambda,
                               family = family, start = start,
                               maxit = maxit, tol = tol)
  }

  return(out)
}

#' Fit a penalised generalised linear model by penalised IRLS on R
#'
#' Fit a generalised linear model with a quadratic penalty on coefficients
#' using an iteratively reweighted least squares (IRLS) approach. The penalty
#' is specified via a diagonal penalty weight vector and a scalar tuning
#' parameter \code{lambda}; the effective penalty matrix is
#' \eqn{\sqrt{\lambda}\,P} where \code{P} = diag(\code{penalty_weights}).
#'
#' @param X **matrix** of predictors with \eqn{n} rows (observations) and
#'   \eqn{p} columns (features).
#' @param y **numeric** response vector of length \eqn{n}.
#' @param weights **numeric** optional observation weights of length \eqn{n}.
#'   Defaults to \code{rep(1, n)}.
#' @param offset **numeric** optional offset vector of length \eqn{n}.
#'   Defaults to \code{rep(0, n)}.
#' @param penalty_weights **numeric** nonnegative vector of length \eqn{p}
#'   giving diagonal entries of the penalty matrix \code{P}. Defaults to
#'   \code{rep(0, p)} (no penalty).
#' @param lambda **numeric** nonnegative scalar that scales the penalty matrix.
#'   Defaults to \code{0} (no penalty).
#' @param family a \pkg{stats}-style family object (for example
#'   \code{gaussian()}, \code{binomial()}, \code{poisson()}). The function
#'   uses \code{family$linkinv}, \code{family$mu.eta} and
#'   \code{family$variance}. Default is \code{gaussian()}.
#' @param start **numeric** optional starting values for the coefficients of
#'   length \eqn{p}. If \code{NULL} (default) coefficients are initialized to
#'   zero. A length mismatch raises an error.
#' @param maxit **integer** maximum number of IRLS iterations. Default is
#'   \code{50}.
#' @param tol **numeric** convergence tolerance on the maximum absolute change
#'   in coefficients between iterations. Default is \code{1e-6}.
#' @param verbose **logical** if \code{TRUE} prints iteration progress.
#'   Default is \code{FALSE}.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{coefficients}}{numeric vector of length \eqn{p} with fitted
#'     regression coefficients.}
#'   \item{\code{eta}}{numeric vector of linear predictors \eqn{X \beta + offset}.}
#'   \item{\code{mu}}{numeric vector of fitted means \eqn{g^{-1}(\eta)} where
#'     \eqn{g^{-1}} is the family link inverse.}
#' }
#'
#' @details
#' The algorithm performs IRLS where, at each iteration, a weighted least
#' squares problem with augmented rows is solved:
#' \preformatted{
#'   A = rbind(sqrt(W) * X, sqrt(lambda) * P)
#'   b = c(sqrt(W) * z, rep(0, p))
#' }
#' where \code{W} are the IRLS weights and \code{z} is the working response.
#' The penalty matrix \code{P} is diagonal with entries given by
#' \code{penalty_weights}. When \code{lambda == 0} or all
#' \code{penalty_weights == 0} the fit reduces to an unpenalised IRLS fit.
#'
#' If \code{start} is provided it must have length equal to the number of
#' columns of \code{X}; otherwise an error is raised.
#'
#' @examples
#' ## Gaussian example
#' set.seed(1)
#' n <- 100; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(1, -0.5, 0, 0, 0)
#' y <- X %*% beta_true + rnorm(n)
#' fit <- penalised_glm_fit(X, y, lambda = 1, penalty_weights = c(0, 1, 1, 1, 1))
#' print(fit$coefficients)
#'
#' ## Binomial example (logit link)
#' set.seed(2)
#' n <- 200; p <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(0.5, -1, 0.2)
#' eta <- X %*% beta_true
#' p_true <- 1 / (1 + exp(-eta))
#' y_bin <- rbinom(n, size = 1, prob = p_true)
#' fit_bin <- penalised_glm_fit_R(X, y_bin, family = binomial(), lambda = 0.5,
#'                              penalty_weights = rep(1, p))
#' head(fit_bin$mu)
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#' @export
penalised_glm_fit_R <- function(X, y, weights = NULL, offset = NULL,
                                penalty_weights = NULL, lambda = 0,
                                family = gaussian(), start = NULL,
                                maxit = 50, tol = 1e-6, verbose = FALSE) {

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset)) offset <- rep(0, n)
  if (is.null(penalty_weights)) penalty_weights <- rep(0, p)
  if (length(penalty_weights) != p) stop("penalty_weights must have length p")

  linkinv  <- family$linkinv
  mu_eta   <- family$mu.eta
  variance <- family$variance

  if (is.null(start)) {
    beta <- rep(0, p)
  } else {
    if (length(start) == p) {
      beta <- as.numeric(start)
    } else {
      stop("incorrect length of starting values")
    }
  }

  eta <- as.vector(X %*% beta + offset)
  mu  <- linkinv(eta)

  converged <- FALSE
  iter_used <- 0

  for (iter in 1:maxit) {
    iter_used <- iter

    mu_eta_val <- mu_eta(eta)
    var_mu <- variance(mu)

    # protect against zero or NA mu_eta or variance
    if (any(!is.finite(mu_eta_val)) || any(!is.finite(var_mu))) {
      warning("non-finite mu_eta or variance encountered; stopping IRLS")
      break
    }

    W <- weights * (mu_eta_val^2 / var_mu)
    # if all weights are zero, break
    if (all(W == 0)) {
      warning("all IRLS weights are zero; returning current coefficients")
      break
    }

    z <- (eta - offset) + (y - mu) / mu_eta_val

    # Weighted normal equations: XtW X and XtW z
    # compute X' (W * X) efficiently
    sqrtW <- sqrt(W)
    WX <- sqrtW * X
    wz <- sqrtW * z

    XtWX <- crossprod(WX)            # p x p
    XtWz <- crossprod(WX, wz)        # p

    # Add penalty: diagonal matrix diag(lambda * penalty_weights)
    if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
      stop("lambda must be a non-negative scalar")
    }
    penalty_diag <- lambda * as.numeric(penalty_weights)
    # NOT DOING THIS: Add small jitter to improve numerical stability when penalty_diag is zero
    #jitter_eps <- 1e-12
    reg_mat <- XtWX
    diag(reg_mat) <- diag(reg_mat) + penalty_diag# + jitter_eps

    # Try Cholesky solve first for speed and stability when PD
    beta_new <- rep(NA_real_, p)
    chol_ok <- FALSE
    chol_err <- NULL
    tryCatch({
      R <- chol(reg_mat)            # R is upper triangular, R' R = reg_mat
      tmp.resp <- backsolve(R, XtWz, transpose = TRUE)  # solves R' tmp.resp = XtWz
      beta_new <- backsolve(R, tmp.resp)                 # solves R beta = tmp.resp
      chol_ok <- TRUE
    }, error = function(e) {
      chol_err <<- conditionMessage(e)
      chol_ok <<- FALSE
    })

    if (!chol_ok) {
      # fallback to QR on the regularized system
      # build small augmented system: [WX; sqrt(penalty_diag) * I] and rhs [wz; 0]
      # but avoid forming huge matrices by using qr on reg_mat directly
      qr_reg <- qr(reg_mat)
      if (qr_reg$rank < p) {
        warning(sprintf("regularized normal matrix rank deficient (rank=%d < p=%d); using qr.coef with tolerance", qr_reg$rank, p))
      }
      beta_new <- qr.coef(qr_reg, XtWz)
      # if qr.coef returns NA, try a damped ridge
      if (any(is.na(beta_new))) {
        damp <- max(1e-8, 1e-6 * sum(diag(reg_mat)))
        diag(reg_mat) <- diag(reg_mat) + damp
        qr_reg2 <- qr(reg_mat)
        beta_new <- qr.coef(qr_reg2, XtWz)
        if (any(is.na(beta_new))) {
          warning("unable to obtain stable solution for penalised GLM; returning previous beta")
          beta_new <- beta
          break
        }
      }
    }

    # check for non-finite solution
    if (any(!is.finite(beta_new))) {
      warning("non-finite coefficients produced; stopping and returning previous coefficients")
      beta_new <- beta
      break
    }

    eta_new <- as.vector(X %*% beta_new + offset)
    mu_new  <- linkinv(eta_new)

    if (max(abs(beta_new - beta), na.rm = TRUE) < tol) {
      beta <- beta_new; eta <- eta_new; mu <- mu_new
      converged <- TRUE
      break
    }

    beta <- beta_new
    eta <- eta_new
    mu  <- mu_new
  }

  list(coefficients = as.numeric(beta),
       # eta = as.numeric(eta),
       # mu = as.numeric(mu),
       converged = converged,
       iterations = iter_used)
}

#' Fit a penalised generalised linear model by penalised IRLS in C++
#'
#' Fit a generalised linear model with a quadratic penalty on coefficients
#' using an iteratively reweighted least squares (IRLS) approach. The penalty
#' is specified via a diagonal penalty weight vector and a scalar tuning
#' parameter \code{lambda}; the effective penalty matrix is
#' \eqn{\sqrt{\lambda}\,P} where \code{P} = diag(\code{penalty_weights}).
#'
#' @param X **matrix** of predictors with \eqn{n} rows (observations) and
#'   \eqn{p} columns (features).
#' @param y **numeric** response vector of length \eqn{n}.
#' @param weights **numeric** optional observation weights of length \eqn{n}.
#'   Defaults to \code{rep(1, n)}.
#' @param offset **numeric** optional offset vector of length \eqn{n}.
#'   Defaults to \code{rep(0, n)}.
#' @param penalty_weights **numeric** nonnegative vector of length \eqn{p}
#'   giving diagonal entries of the penalty matrix \code{P}. Defaults to
#'   \code{rep(0, p)} (no penalty).
#' @param lambda **numeric** nonnegative scalar that scales the penalty matrix.
#'   Defaults to \code{0} (no penalty).
#' @param family a \pkg{stats}-style family object (for example
#'   \code{gaussian()}, \code{binomial()}, \code{poisson()}). The function
#'   uses \code{family$linkinv}, \code{family$mu.eta} and
#'   \code{family$variance}. Default is \code{gaussian()}.
#' @param start **numeric** optional starting values for the coefficients of
#'   length \eqn{p}. If \code{NULL} (default) coefficients are initialised to
#'   zero. A length mismatch raises an error.
#' @param maxit **integer** maximum number of IRLS iterations. Default is
#'   \code{50}.
#' @param tol **numeric** convergence tolerance on the maximum absolute change
#'   in coefficients between iterations. Default is \code{1e-6}.
#' @param verbose **logical** if \code{TRUE} prints iteration progress.
#'   Default is \code{FALSE}.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{coefficients}}{numeric vector of length \eqn{p} with fitted
#'     regression coefficients.}
#'   \item{\code{eta}}{numeric vector of linear predictors \eqn{X \beta + offset}.}
#'   \item{\code{mu}}{numeric vector of fitted means \eqn{g^{-1}(\eta)} where
#'     \eqn{g^{-1}} is the family link inverse.}
#' }
#'
#' @details
#' The algorithm performs IRLS where, at each iteration, a weighted least
#' squares problem with augmented rows is solved:
#' \preformatted{
#'   A = rbind(sqrt(W) * X, sqrt(lambda) * P)
#'   b = c(sqrt(W) * z, rep(0, p))
#' }
#' where \code{W} are the IRLS weights and \code{z} is the working response.
#' The penalty matrix \code{P} is diagonal with entries given by
#' \code{penalty_weights}. When \code{lambda == 0} or all
#' \code{penalty_weights == 0} the fit reduces to an unpenalised IRLS fit.
#'
#' If \code{start} is provided it must have length equal to the number of
#' columns of \code{X}; otherwise an error is raised.
#'
#' @examples
#' ## Gaussian example
#' set.seed(1)
#' n <- 100; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(1, -0.5, 0, 0, 0)
#' y <- X %*% beta_true + rnorm(n)
#' fit <- penalised_glm_fit(X, y, lambda = 1, penalty_weights = c(0, 1, 1, 1, 1))
#' print(fit$coefficients)
#'
#' ## Binomial example (logit link)
#' set.seed(2)
#' n <- 200; p <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' beta_true <- c(0.5, -1, 0.2)
#' eta <- X %*% beta_true
#' p_true <- 1 / (1 + exp(-eta))
#' y_bin <- rbinom(n, size = 1, prob = p_true)
#' fit_bin <- penalised_glm_fit(X, y_bin, family = binomial(), lambda = 0.5,
#'                              penalty_weights = rep(1, p))
#' head(fit_bin$mu)
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#' @export
penalised_glm_fit_C <- function(X, y, weights = NULL, offset = NULL,
                              penalty_weights = NULL, lambda = 0,
                              family = gaussian(), start = NULL,
                              maxit = 50, tol = 1e-6
) {

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset)) offset <- rep(0, n)
  if (is.null(penalty_weights)) penalty_weights <- rep(0, p)
  if (length(penalty_weights) != p) stop("penalty_weights must have length p")

  if (is.null(start)) {
    start <- rep(0, p)
  } else {
    if (length(start) == p) {
      start <- as.numeric(start)
    } else {
      stop("incorrect length of starting values")
    }
  }

  family = switch(family$family,
                  poisson = 0,
                  binomial = 1,
                  gaussian = 2)

  res <- .Call(
    "penalised_glm_irls_cpp",
    X, y, weights, offset,
    penalty_weights, lambda, start,
    family, maxit, tol
  )

  list(
    coefficients = res$coefficients,
    converged = res$converged,
    iterations = res$iterations
  )
}

#' E-step: posterior probabilities for Correlated Species Archetype Models (CSAM)
#'
#' Compute posterior membership probabilities for each species to archetypes
#' in the E-step of the EM algorithm used by the \pkg{csam} package. For each
#' species \(j\) and archetype \(k\) the function evaluates the archetype-
#' specific GLM log-likelihood under the supplied \code{family}, adds the log
#' mixing proportion \code{log(pi[k])}, and normalizes across archetypes to
#' produce posterior probabilities \eqn{\tau_{jk}}. An observation-level
#' offset is constructed from latent site covariates \code{U} and species
#' loadings \code{Lambda} to capture residual correlation among species.
#'
#' @param Y numeric matrix of species responses with \eqn{n} rows (sites) and
#'   \eqn{s} columns (species).
#' @param X numeric design matrix with \eqn{n} rows (sites) and \eqn{p}
#'   columns (environmental covariates).
#' @param beta0 numeric vector of length \eqn{s} of species-specific intercepts.
#' @param B numeric matrix of archetype-specific coefficients with \eqn{g}
#'   rows and \eqn{p} columns (row \code{k} = coefficients for archetype \code{k}).
#' @param U numeric matrix of latent site covariates with \eqn{n} rows and
#'   \eqn{r} columns.
#' @param Lambda numeric matrix of species loadings with \eqn{s} rows and
#'   \eqn{r} columns; row \code{j} contains loadings for species \code{j}.
#' @param phi numeric dispersion/scale parameter (kept for API compatibility;
#'   not used directly by this function).
#' @param pi numeric vector of length \eqn{g} of prior archetype mixing
#'   proportions (positive, need not be normalized).
#' @param family a \pkg{stats}-style family object (e.g., \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}). The function uses
#'   \code{family$linkinv}, \code{family$dev.resids} and \code{family$aic}.
#' @param trunc.interval logical; if \code{TRUE} apply a deterministic
#'   truncation/rescaling transform to posterior probabilities (default
#'   \code{FALSE}).
#'
#' @return A numeric matrix of dimension \eqn{s \times g} (species by
#'   archetypes). Row \code{j} contains posterior probabilities
#'   \eqn{\tau_{j1},\dots,\tau_{jg}} for species \code{j}. Rows sum to 1
#'   unless \code{trunc.interval = TRUE}, in which case a deterministic
#'   rescaling is applied columnwise.
#'
#' @details
#' In the CSAM notation used by \pkg{csam}, the linear predictor for species
#' \code{j} under archetype \code{k} is
#' \deqn{\eta_{jk} = \beta_{0j} + X B_{k} + U \lambda_{j},}
#' where \eqn{\lambda_j} is row \code{j} of \code{Lambda}. Per-observation
#' deviance/AIC-based log-likelihood contributions are obtained from
#' \code{family} and summed across sites to form \code{loglik_mat[j,k]}. The
#' posterior log-probabilities are \code{loglik_mat + log(pi)} and are
#' stabilized by subtracting the row-wise maximum before exponentiating and
#' normalizing.
#'
#' The optional \code{trunc.interval} transform rescales each column of the
#' posterior matrix using the fixed formula implemented in the code; this is
#' a heuristic post-processing step used in some CSAM fitting routines.
#'
#' @examples
#' set.seed(1)
#' n <- 60; s <- 5; p <- 3; r <- 2; g <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' U <- matrix(rnorm(n * r), n, r)
#' B <- matrix(rnorm(g * p), g, p)
#' beta0 <- rnorm(s)
#' Lambda <- matrix(rnorm(s * r), s, r)
#' pi <- c(0.7, 0.3)
#' Y <- matrix(rnorm(n * s), n, s)
#' tau <- estep_post_probs(Y, X, beta0, B, U, Lambda, phi = 1, pi,
#'                         family = gaussian(), trunc.interval = FALSE)
#' dim(tau)  # s x g: species by archetypes
#'
#' @seealso \code{\link[stats]{family}}
#' @export
estep_post_probs <- function(Y, X, par.list, family,
                                  trunc.interval = FALSE) {

  beta0  <- par.list$beta0
  B      <- par.list$B
  pi     <- par.list$pi
  U      <- par.list$U
  Lambda <- par.list$Lambda

  n <- nrow(Y)
  s <- ncol(Y)
  g <- nrow(B)

  ## ---- Shared linear predictor components (computed once) ----
  XB <- X %*% t(B)          # n × g
  UL <- U %*% t(Lambda)     # n × s

  loglik_mat <- matrix(0, s, g)

  ## ---- E-step loops (species-level, parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    eta_base <- beta0[j] + UL[, j]

    for (k in seq_len(g)) {

      eta <- eta_base + XB[, k]
      mu  <- family$linkinv(eta)

      ## deviances already include all family-specific terms
      dev <- family$dev.resids(yj, mu, wt = rep(1, n))

      ## relative log-likelihood (constants drop out)
      loglik_mat[j, k] <- -0.5 * sum(dev)
    }
  }

  ## ---- Posterior normalisation (softmax) ----
  logpost <- sweep(loglik_mat, 2, log(pi), "+")

  ## stabilise row-wise
  row_max <- apply(logpost, 1, max)
  post_unnorm <- exp(logpost - row_max)
  tau <- post_unnorm / rowSums(post_unnorm)

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

#' M-step: update archetype (archetype-level) regression coefficients for CSAM
#'
#' Update archetype-specific regression coefficients \(B\) in the conditional
#' M-step of the Correlated Species Archetype Model (CSAM) EM algorithm.
#' Given current species intercepts, site latent variables and species
#' loadings, and posterior species-to-archetype probabilities \eqn{\tau},
#' this function fits a weighted GLM for each archetype by stacking species
#' observations and using \eqn{\tau_{jk}} as observation weights.
#'
#' @param Y numeric matrix of species responses with \eqn{n} rows (sites) and
#'   \eqn{s} columns (species).
#' @param X numeric design matrix with \eqn{n} rows (sites) and \eqn{p}
#'   columns (environmental covariates).
#' @param beta0 numeric vector of length \eqn{s} of species-specific intercepts.
#' @param B numeric matrix of current archetype coefficients with \eqn{g}
#'   rows and \eqn{p} columns; row \code{k} contains coefficients for archetype \code{k}.
#' @param U numeric matrix of latent site covariates with \eqn{n} rows and
#'   \eqn{r} columns.
#' @param Lambda numeric matrix of species loadings with \eqn{s} rows and
#'   \eqn{r} columns.
#' @param phi numeric dispersion/scale parameter (kept for API compatibility;
#'   not used directly by this function).
#' @param pi numeric vector of length \eqn{g} of prior archetype mixing
#'   proportions (not used directly by this update but part of the model state).
#' @param tau numeric matrix of posterior species-to-archetype probabilities
#'   with dimension \eqn{s \times g} (rows = species, columns = archetypes).
#' @param family a \pkg{stats}-style family object (e.g., \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}) used for the GLM fits.
#' @param maxit integer maximum number of iterations passed to the internal
#'   \code{\link[stats]{glm.fit}} call (default \code{1}). A small number of
#'   iterations is typical because this is a conditional update inside EM.
#'
#' @return A numeric matrix with the same dimensions as \code{B} containing
#'   the updated archetype regression coefficients.
#'
#' @details
#' The update proceeds archetype-by-archetype. For archetype \code{k} the
#' function stacks all species' responses and design matrices into a single
#' long vector/matrix of length \eqn{n s} and fits a weighted GLM with
#' observation-level offsets equal to \eqn{\beta_{0j} + U \lambda_j} for the
#' block of rows corresponding to species \code{j}. The weights for rows
#' corresponding to species \code{j} are given by \eqn{\tau_{jk}}. The fitted
#' coefficients for archetype \code{k} replace the corresponding row of
#' \code{B}.
#'
#' The implementation uses \code{stats::glm.fit} and suppresses warnings
#' about non-convergence because only a small number of iterations is
#' typically requested inside the EM loop.
#'
#' @examples
#' set.seed(1)
#' n <- 50; s <- 4; p <- 3; r <- 2; g <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' U <- matrix(rnorm(n * r), n, r)
#' B <- matrix(rnorm(g * p), g, p)
#' beta0 <- rnorm(s)
#' Lambda <- matrix(rnorm(s * r), s, r)
#' pi <- c(0.6, 0.4)
#' Y <- matrix(rnorm(n * s), n, s)
#' tau <- matrix(runif(s * g), s, g)
#' tau <- tau / rowSums(tau)
#' B_upd <- mstep_arch_pars(Y, X, beta0, B, U, Lambda, phi = 1, pi, tau,
#'                          family = gaussian(), maxit = 1)
#'
#' @seealso \code{\link[stats]{glm.fit}}, \code{\link{estep_post_probs}}
#' @export
mstep_arch_pars <- function(Y, X, par.list, tau,
                            family, maxit = 1) {

  beta0 = par.list$beta0; B = par.list$B; pi = par.list$pi; U = par.list$U; Lambda = par.list$Lambda; phi = par.list$phi

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X)
  B_new <- B

  # note the following removes non-convergence as we are typically using only a single iteration
  safe_glm_fit <- function(...,
                           pattern = "algorithm did not converge") {
    withCallingHandlers(
      stats::glm.fit(...),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl(pattern, msg, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        } else {
          # re-emit other warnings
          warning(w)
        }
      }
    )
  }

  ## ---- Shared linear predictor (fixed in this CM-step) ----
  UL <- U %*% t(Lambda)    # n × s

  ## ---- Pre-build stacked design matrix ONCE ----
  # Stack X by species
  X_stack <- X[rep(seq_len(n), times = s), ]

  ns <- n * s
  y_stack       <- numeric(ns)
  offset_stack  <- numeric(ns)
  weights_stack <- numeric(ns)

  for (k in seq_len(g)) {

    ## ---- Fill stacked vectors by species ----
    for (j in seq_len(s)) {
      idx <- ((j - 1) * n + 1):(j * n)

      y_stack[idx]       <- Y[, j]
      offset_stack[idx]  <- beta0[j] + UL[, j]
      weights_stack[idx] <- tau[j, k]
    }

    fit <- safe_glm_fit(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    )

    if (is.null(fit$coefficients) || any(is.na(fit$coefficients))) {
      warning(sprintf("glm.fit failed for archetype %d; leaving B[k,] unchanged", k))
    } else {
      B_new[k, ] <- fit$coefficients
    }

    # fit <- suppressWarnings(stats::glm.fit(
    #   x = X_stack,
    #   y = y_stack,
    #   weights = weights_stack,
    #   offset = offset_stack,
    #   family = family,
    #   control = list(maxit = maxit)
    # ))

    # B_new[k, ] <- fit$coefficients
  }

  B_new
}

#' M-step: update species-level parameters for CSAM (species intercepts, loadings, dispersion)
#'
#' Perform the conditional M-step that updates species-specific parameters in
#' the Correlated Species Archetype Model (CSAM). For each species \(j\) the
#' function stacks observations across archetypes and fits a weighted,
#' penalised GLM to update the species intercept \eqn{\beta_{0j}} and species
#' loadings \eqn{\lambda_j} (row \code{j} of \code{Lambda}). The posterior
#' species-to-archetype probabilities \eqn{\tau_{jk}} are used as observation
#' weights and the archetype fixed-effects contribution is included as an
#' offset.
#'
#' @param Y numeric matrix of species responses with \eqn{n} rows (sites) and
#'   \eqn{s} columns (species).
#' @param X numeric design matrix with \eqn{n} rows (sites) and \eqn{p}
#'   columns (environmental covariates). Used to compute archetype offsets.
#' @param B numeric matrix of archetype regression coefficients with \eqn{g}
#'   rows and \eqn{p} columns (row \code{k} = coefficients for archetype \code{k}).
#' @param U numeric matrix of latent site covariates with \eqn{n} rows and
#'   \eqn{d} columns (site-level latent variables).
#' @param beta0 numeric vector of current species intercepts of length \eqn{s}.
#' @param Lambda numeric matrix of current species loadings with \eqn{s} rows
#'   and \eqn{d} columns (row \code{j} = loadings for species \code{j}).
#' @param phi numeric vector or scalar of current dispersion/scale parameters
#'   (one per species when applicable).
#' @param pi numeric vector of archetype mixing proportions (length \eqn{g});
#'   included for API consistency but not used directly in this update.
#' @param tau numeric matrix of posterior species-to-archetype probabilities
#'   with dimension \eqn{s \times g} (rows = species, columns = archetypes).
#' @param family a \pkg{stats}-style family object (e.g., \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}) used for the GLM fits.
#' @param psi2 numeric nonnegative scalar penalty tuning parameter applied to
#'   species loadings (ridge penalty). Default \code{0} (no penalty).
#' @param maxit integer maximum number of IRLS iterations passed to the
#'   internal \code{\link{penalised_glm_fit}} call (default \code{1}).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{beta0}}{numeric vector of updated species intercepts (length \eqn{s}).}
#'   \item{\code{Lambda}}{numeric matrix of updated species loadings (dimension \eqn{s \times d}).}
#'   \item{\code{phi}}{numeric vector of updated dispersion/scale parameters (length \eqn{s}) when applicable.}
#' }
#'
#' @details
#' For each species \code{j} the function constructs a stacked dataset of
#' length \eqn{n g} by repeating the site responses and design for each
#' archetype. The stacked design has an intercept column and the site latent
#' variables \code{U} as covariates; the archetype contribution
#' \eqn{X B_k}, is included as an observation-level offset for rows
#' corresponding to archetype \code{k}. Observation weights are given by
#' \eqn{\tau_{jk}}. A penalised GLM is fitted via \code{\link{penalised_glm_fit}}
#' with penalty weights \code{c(0, rep(1, d))} so that the intercept is
#' unpenalised while the loadings receive a ridge penalty scaled by
#' \code{psi2}. When \code{family$family == "gaussian"} the species-specific
#' dispersion \eqn{\phi_j} is re-estimated as the weighted mean squared error
#' using the IRLS weights.
#'
#' The function returns the updated \code{beta0}, \code{Lambda} and \code{phi}
#' which can be used in subsequent EM iterations.
#'
#' @examples
#' set.seed(42)
#' n <- 40; s <- 6; p <- 2; d <- 2; g <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' U <- matrix(rnorm(n * d), n, d)
#' B <- matrix(rnorm(g * p), g, p)
#' beta0 <- rnorm(s)
#' Lambda <- matrix(rnorm(s * d), s, d)
#' pi <- rep(1/g, g)
#' Y <- matrix(rnorm(n * s), n, s)
#' tau <- matrix(runif(s * g), s, g); tau <- tau / rowSums(tau)
#' res <- mstep_species_pars(Y, X, B, U, beta0, Lambda, phi = rep(1, s),
#'                           pi, tau, family = gaussian(), psi2 = 0.1, maxit = 2)
#' str(res)
#'
#' @seealso \code{\link{penalised_glm_fit}}, \code{\link{mstep_arch_pars}},
#'   \code{\link{estep_post_probs}}
#' @export
mstep_species_pars <- function(
    Y, X, par.list, tau,
    family, psi2 = 0, maxit = 1, backend = c("C++", "R")
) {
  backend <- match.arg(backend)

  beta0  <- par.list$beta0
  B      <- par.list$B
  U      <- par.list$U
  Lambda <- par.list$Lambda
  phi    <- par.list$phi

  n <- nrow(Y)
  s <- ncol(Y)
  g <- nrow(B)
  d <- ncol(U)

  ## ---- Shared quantities (once per ECM iteration) ----
  XB <- X %*% t(B)          # n × g
  XU <- cbind(1, U)         # n × (d+1)

  ## ---- Allocate stacking buffers ONCE ----
  ng <- n * g
  y_stack       <- numeric(ng)
  offset_stack  <- numeric(ng)
  weights_stack <- numeric(ng)
  X_stack       <- matrix(0, nrow = ng, ncol = d + 1)

  beta0_new  <- beta0
  Lambda_new <- Lambda
  phi_new    <- phi

  ## ---- Species loop (parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]

    for (k in seq_len(g)) {
      idx <- ((k - 1) * n + 1):(k * n)

      y_stack[idx]       <- yj
      X_stack[idx, ]     <- XU
      offset_stack[idx]  <- XB[, k]
      weights_stack[idx] <- tau[j, k]
    }

    fit <- penalised_glm_fit(
      X       = X_stack,
      y       = y_stack,
      weights = weights_stack,
      offset  = offset_stack,
      penalty_weights = c(0, rep(1, d)),
      lambda  = psi2,
      family  = family,
      start   = c(beta0_new[j], Lambda_new[j, ]),
      maxit   = maxit,
      backend = backend
    )

    if (all(is.finite(fit$coefficients))) {
      beta0_new[j]    <- fit$coefficients[1]
      Lambda_new[j, ] <- fit$coefficients[-1]
    }

    if (family$family == "gaussian") {
      phi_new[j] <- sum(weights_stack * (y_stack - fit$mu)^2) /
        sum(weights_stack)
    }
  }

  list(beta0 = beta0_new, Lambda = Lambda_new, phi = phi_new)
}


#' M-step: update site scores (latent site covariates) for CSAM
#'
#' Update site-level latent variables \code{U} (site scores) in the conditional
#' M-step of the Correlated Species Archetype Model (CSAM). For each site the
#' function stacks species-by-archetype observations and fits a penalised GLM
#' to update the site's latent coordinates. Species loadings \code{Lambda}
#' serve as covariates in the per-site regression, archetype contributions
#' \code{X %*% B[k, ]} enter as offsets, and posterior species-to-archetype
#' probabilities \eqn{\tau_{jk}} provide observation weights.
#'
#' @param Y **numeric matrix** of species responses with \eqn{n} rows (sites)
#'   and \eqn{s} columns (species).
#' @param X **numeric matrix** design matrix with \eqn{n} rows (sites) and
#'   \eqn{p} columns (environmental covariates).
#' @param B **numeric matrix** of archetype regression coefficients with
#'   \eqn{g} rows and \eqn{p} columns (row \code{k} = coefficients for archetype \code{k}).
#' @param beta0 **numeric** vector of species intercepts of length \eqn{s}.
#' @param Lambda **numeric matrix** of species loadings with \eqn{s} rows and
#'   \eqn{d} columns (row \code{j} = loadings for species \code{j}).
#' @param U **numeric matrix** of current latent site covariates with
#'   \eqn{n} rows and \eqn{d} columns (row \code{i} = site \code{i} scores).
#' @param phi **numeric** dispersion/scale parameter (kept for API compatibility;
#'   not used directly by this function).
#' @param pi **numeric** vector of archetype mixing proportions (length \eqn{g});
#'   included for API consistency but not used directly in this update.
#' @param tau **numeric matrix** of posterior species-to-archetype probabilities
#'   with dimension \eqn{s \times g} (rows = species, columns = archetypes).
#' @param family a \pkg{stats}-style family object (for example
#'   \code{gaussian()}, \code{binomial()}, \code{poisson()}). The function
#'   uses \code{family$linkinv}, \code{family$dev.resids} and related methods
#'   via the internal GLM solver.
#' @param psi1 **numeric** nonnegative scalar penalty tuning parameter applied
#'   to site scores (ridge penalty). Default \code{0} (no penalty).
#' @param maxit **integer** maximum number of IRLS iterations passed to the
#'   internal \code{\link{penalised_glm_fit}} call (default \code{1}).
#'
#' @return A numeric matrix of the same dimension as \code{U} (\eqn{n \times d})
#'   containing the updated site scores (row \code{i} = updated scores for site \code{i}).
#'
#' @details
#' For each site \code{i} the function constructs a stacked dataset of length
#' \eqn{s g} by iterating over species \code{j} and archetypes \code{k}. The
#' stacked design matrix uses the species loadings \eqn{\lambda_j} (row
#' \code{j} of \code{Lambda}) as covariates; the archetype contribution
#' \eqn{X[i, ] B_k}, and species intercept \eqn{\beta_{0j}} are included as
#' an observation-level offset. Observation weights are given by
#' \eqn{\tau_{jk}}. A penalised GLM is fitted via \code{\link{penalised_glm_fit}}
#' with penalty weights \code{rep(1, d)} so that all dimensions of the site
#' score receive a ridge penalty scaled by \code{psi1}. The fitted
#' coefficients replace the corresponding row of \code{U}.
#'
#' This conditional update treats each site independently and is typically
#' used inside an outer EM loop that alternates E- and M-steps for CSAM.
#'
#' @examples
#' set.seed(123)
#' n <- 30; s <- 5; p <- 2; d <- 2; g <- 3
#' X <- matrix(rnorm(n * p), n, p)
#' B <- matrix(rnorm(g * p), g, p)
#' beta0 <- rnorm(s)
#' Lambda <- matrix(rnorm(s * d), s, d)
#' U <- matrix(rnorm(n * d), n, d)
#' pi <- rep(1/g, g)
#' Y <- matrix(rnorm(n * s), n, s)
#' tau <- matrix(runif(s * g), s, g); tau <- tau / rowSums(tau)
#' U_upd <- mstep_site_scores(Y, X, B, beta0, Lambda, U, phi = 1, pi, tau,
#'                            family = gaussian(), psi1 = 0.1, maxit = 2)
#' dim(U_upd)  # n x d
#'
#' @seealso \code{\link{penalised_glm_fit}}, \code{\link{mstep_arch_pars}},
#'   \code{\link{mstep_species_pars}}, \code{\link{estep_post_probs}}
#' @keywords csam site mstep latent penalised
#' @export
mstep_site_scores <- function(Y, X, par.list, tau,
                              family, psi1 = 0, maxit = 1, backend = c("C++", "R")) {
  backend <- match.arg(backend)

  beta0  <- par.list$beta0
  B      <- par.list$B
  U      <- par.list$U
  Lambda <- par.list$Lambda

  n <- nrow(Y)
  s <- ncol(Y)
  g <- nrow(B)
  d <- ncol(Lambda)

  ## ---- Shared linear predictor (fixed in this CM-step) ----
  XB <- X %*% t(B)   # n × g

  U_new <- U

  sg <- s * g

  ## ---- Allocate stacks ONCE ----
  y_stack       <- numeric(sg)
  offset_stack  <- numeric(sg)
  weights_stack <- numeric(sg)
  X_stack       <- matrix(0, nrow = sg, ncol = d)

  penalty_weights <- rep(1, d)

  for (i in seq_len(n)) {

    row <- 1L

    ## ---- Fill stacks by species × archetype ----
    for (j in seq_len(s)) {
      lj <- Lambda[j, ]
      bj <- beta0[j]

      for (k in seq_len(g)) {
        y_stack[row]       <- Y[i, j]
        X_stack[row, ]     <- lj
        offset_stack[row]  <- bj + XB[i, k]
        weights_stack[row] <- tau[j, k]
        row <- row + 1L
      }
    }

    fit <- penalised_glm_fit(
      X = X_stack,
      y = y_stack,
      weights = weights_stack,
      offset = offset_stack,
      penalty_weights = penalty_weights,
      lambda = psi1,
      family = family,
      start = U_new[i, ],
      maxit = maxit,
      backend = backend
    )

    if (!is.null(fit$coefficients) &&
        all(is.finite(fit$coefficients))) {
      U_new[i, ] <- fit$coefficients
    } else {
      warning(sprintf(
        "penalised_glm_fit failed for site %d; leaving U[i, ] unchanged",
        i
      ))
    }
  }

  U_new
}


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

      dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      loglik_mat[j, k] <- -0.5 * sum(dev)
    }
  }

  ## ---- Posterior normalisation ----
  logpost <- sweep(loglik_mat, 2, log(pi), "+")
  row_max <- apply(logpost, 1, max)
  post_unnorm <- exp(logpost - row_max)
  tau <- post_unnorm / rowSums(post_unnorm)

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
                            family, maxit = 1) {

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
mstep0_species_pars <- function(Y, X, par.list, tau,
                               family, maxit = 1) {

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
  if (is.null(object$U)) {
    U <- matrix(rep(0, nrow(Y)), ncol = 1)
  } else {
    U <- object$U
  }
  if(is.null(object$Lambda)) {
    Lambda <- matrix(rep(0, ncol(Y)), ncol = 1)
  } else {
    Lambda <- object$Lambda
  }
  phi <- object$phi
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

init.fa.pars_gllvm <- function(Y, X, g = 3, family = poisson(), d = 2) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise usual SAM parameters
  start.pars = csam::init.sam.pars(Y, X, g = g, family = family)

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013
  res = matrix(rep(0, n * s), n, s)
  for (sp in 1:s) {
    if (family$family == "binomial") {
      tmp.offset = start.pars$beta0[sp] + X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
      tmp.m = stats::glm(y ~ 0, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
      res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m), quantileFunction = qnorm)
    } else {
      tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
      tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
      res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
    }
  }
  # perform factor analysis on the residuals
  tmp.fa = gllvm::gllvm(y = res, family = gaussian(), num.lv = d)
  tmp.dat = data.frame(y = as.vector(res), x = matrix(rep(1, n * s), ncol = 1))

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
#' @importFrom DHARMa simulateResiduals
#' @importFrom stats glm
init.fa.pars <- function(Y, X, g = 3, family = poisson(), d = 2) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise usual SAM parameters
  start.pars = csam::init.sam.pars(Y, X, g = g, family = family)

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013

  # with hard class labelling:
  res = matrix(rep(0, n * s), n, s)
  for (sp in 1:s) {
    if (family$family == "binomial") {
      tmp.offset = start.pars$beta0[sp] + X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
      tmp.m = stats::glm(y ~ 0, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
      res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m), quantileFunction = qnorm)
    } else {
      tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
      tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
      res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
    }
    # res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth") # CON: uses an unexported function, throws error for binomial, requires glmmTMB() rather than glm() above
    # res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m), quantileFunction = qnorm) # CON: uses a non-exported function
    # res[ , sp] = DHARMa::simulateResiduals(tmp.m)$scaledResiduals
    # res[ , sp] = statmod::qresiduals(tmp.m) # CON: leads to Inf/-Inf values
  }

  # # with weights for each archetype:
  # res = matrix(rep(0, n * s * g), n * g, s)
  # for (j in 1:s) {
  #   y_stack <- rep(NA, n * g)
  #   weights_stack <- rep(0, n * g)
  #   offset_stack <- rep(0, n * g)
  #   for (k in 1:g) {
  #     idx <- ((k - 1) * n + 1):(k * n)
  #     y_stack[idx] <- Y[, j]
  #     offset_stack[idx] <- as.vector(X %*% start.pars$B[k, ])
  #     weights_stack[idx] <- start.pars$pi[k]
  #   }
  #   # tmp.m = stats::glm(y ~ 0, data = data.frame(y = y_stack), family = family, offset = offset_stack, weights = weights_stack) # CON: see options for qresiduals below
  #   if (family$family == "binomial") {
  #     tmp.m = glmmTMB::glmmTMB(cbind(y, 1 - y) ~ 1, data = data.frame(y = y_stack), family = family, offset = offset_stack, weights = weights_stack)
  #     res[ , j] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m), quantileFunction = qnorm)
  #   } else {
  #     tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = y_stack), family = family, offset = offset_stack, weights = weights_stack)
  #     res[ , j] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
  #   }
  # }

  # put data in long format
  tmp.dat = data.frame(res = as.vector(res), sp = factor(rep(1:s, each = nrow(res))), site = factor(rep(1:n, s)))
  # tmp.dat = data.frame(res = as.vector(res), sp = factor(rep(1:s, each = nrow(res))), site = factor(rep(1:n, s * g)), wt = rep(rep(start.pars$pi, each = n), s))
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

#' Calculate the log-likelihood marginalised over the true archetype labels (internal)
#'
#' @keywords internal
mzll <- function(Y, X, par.list, g = 3, d = 2, family = poisson(),
                 psi1 = 0, psi2 = 0) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  beta0 <- if (!is.null(par.list$beta0)) par.list$beta0 else rep(0, s)
  B      <- if (!is.null(par.list$B))     par.list$B     else matrix(rep(0, g * p), g, p)
  pi     <- if (!is.null(par.list$pi))    par.list$pi    else rep(1 / g, g)
  U      <- if (!is.null(par.list$U))     par.list$U     else matrix(rep(0, n * d), n, d)
  Lambda <- if (!is.null(par.list$Lambda)) par.list$Lambda else matrix(rep(0, s * d), s, d)
  phi    <- if (!is.null(par.list$phi))   par.list$phi   else rep(1, s)

  # some checks on warm starts
  if (g != nrow(B)) {
    stop("mismatch between archetype parameters, B, and specified g")
  }
  if (p != ncol(B)) {
    stop("mismatch between predictor parameters, B, and columns in supplied X")
  }
  if (s != length(beta0)) {
    stop("mismatch between response intercepts, beta0, and response columns in supplied Y")
  }
  if (s != length(phi)) {
    stop("mismatch between response dispersion parameters, phi, and response columns in supplied Y")
  }
  if (g != length(pi)) {
    stop("mismatch between mixing parameters, pi, and specified g")
  }
  if (s != nrow(Lambda)) {
    stop("mismatch between factor loadings, Lambda, and response columns in supplied Y")
  }
  if (n != nrow(U)) {
    stop("mismatch between factor scores, U, and response rows in supplied Y")
  }
  if (d != ncol(Lambda)) {
    stop("mismatch between factor loadings, Lambda, and specified d")
  }
  if (d != ncol(U)) {
    stop("mismatch between factor scores, U, and specified d")
  }

  ## ---- Shared linear predictor components ----
  XB <- X %*% t(B)          # n × g
  UL <- U %*% t(Lambda)     # n × s

  ll <- 0

  ## ---- Species loop (parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    eta_base <- beta0[j] + UL[, j]
    comp <- numeric(g)

    for (k in seq_len(g)) {

      eta <- eta_base + XB[, k]
      mu  <- family$linkinv(eta)

      dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      comp[k] <- -0.5 * sum(dev) + log(pi[k])
    }

    ## numerically stable log-sum-exp
    m <- max(comp)
    ll <- ll + m + log(sum(exp(comp - m)))
  }

  ## ---- Penalised log-likelihood ----
  pll <- ll -
    0.5 * psi1 * sum(U^2) -
    0.5 * psi2 * sum(Lambda^2)

  pll
}


#' Calculate the fitted values of a CSAM/SAM, \eqn{\mu_{ijk}}, and compute the posterior predictive mean (internal)
#'
#' @keywords internal
fitted_csam_internal <- function(par.list, X, tau, family,
                           type = c("response", "link"),
                           archetype_specific = FALSE) {

  type <- match.arg(type)

  beta0  <- par.list$beta0      # length s
  B      <- par.list$B          # g × p
  U      <- par.list$U          # n × d or NULL
  Lambda <- par.list$Lambda     # s × d or NULL

  n <- nrow(X)
  s <- length(beta0)
  g <- nrow(B)

  # --- linear predictor components ---
  XB <- X %*% t(B)   # n × g

  # factor-analytic term (correlated CSAM)
  if (!is.null(U) && !is.null(Lambda)) {
    UL <- U %*% t(Lambda)   # n × s
  } else {
    UL <- matrix(0, n, s)   # standard SAM
  }

  # allocate eta: n × s × g
  eta <- array(0, dim = c(n, s, g))

  # fill eta for each archetype
  for (k in 1:g) {
    # eta_{ijk} = beta0_j + XB_{ik} + UL_{ij}
    eta[,,k] <- sweep(UL, 2, beta0, "+") + XB[,k]
  }

  # return archetype-specific linear predictors if requested
  if (archetype_specific && type == "link") {
    return(eta)
  }

  # inverse link
  invlink <- family$linkinv
  mu <- invlink(eta)   # n × s × g

  # return archetype-specific responses if requested
  if (archetype_specific && type == "response") {
    return(mu)
  }

  # --- posterior predictive mean (default) ---
  # response scale:
  #   fitted_{ij} = sum_k tau[j,k] * mu_{ijk}
  #
  # link scale:
  #   fitted_{ij} = sum_k tau[j,k] * eta_{ijk}

  if (type == "response") {
    out <- matrix(0, n, s)
    for (k in 1:g) {
      out <- out + mu[,,k] * rep(tau[,k], each = n)
    }
    return(out)
  }

  if (type == "link") {
    out <- matrix(0, n, s)
    for (k in 1:g) {
      out <- out + eta[,,k] * rep(tau[,k], each = n)
    }
    return(out)
  }
}

#' Calculate the posterior membership probabilities for a model fitted via TMB (internal)
#'
#' @keywords internal
tau_from_tmb_fit <- function(Y, X, par.list, family) {
  n <- nrow(Y); s <- ncol(Y)
  beta0  <- par.list$beta0
  B      <- par.list$B      # g × p
  U      <- par.list$U
  pi     <- par.list$pi     # length g
  Lambda <- par.list$Lambda # s × d
  phi    <- par.list$phi
  g      <- nrow(B)

  loglik <- matrix(0, s, g)

  # --- linear predictor components ---
  XB <- X %*% t(B)   # n × g

  # factor-analytic term (correlated CSAM)
  if (!is.null(U) && !is.null(Lambda)) {
    UL <- U %*% t(Lambda)   # n × s
  } else {
    UL <- matrix(0, n, s)   # standard SAM
  }

  ## ---- Species-level E-step (parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    eta_base <- beta0[j] + UL[, j]

    for (k in seq_len(g)) {

      eta <- eta_base + XB[, k]
      mu  <- family$linkinv(eta)

      dev <- family$dev.resids(yj, mu, wt = rep(1, n))
      loglik[j, k] <- -0.5 * sum(dev)
    }
  }

  ## ---- Add log pi and normalise ----
  logpost <- sweep(loglik, 2, log(pi), "+")
  row_max <- apply(logpost, 1, max)
  post_unnorm <- exp(logpost - row_max)
  tau <- post_unnorm / rowSums(post_unnorm)

  tau
}
