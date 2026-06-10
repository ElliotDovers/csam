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
                                maxit = 50, tol = 1e-8, verbose = FALSE) {

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
                              maxit = 50, tol = 1e-8
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
    iterations = res$iterations,
    boundary = res$boundary
  )
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
init.sam.pars <- function(Y, X, g = 3, family = poisson(), update.intercepts = FALSE, mixing.props.equal = FALSE) {

  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise parameter list
  start.pars = list(beta0 = rep(0, s), B = matrix(rep(0, g * p), g, p), pi = rep(1/g, g), phi = rep(1, s))

  # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013
  beta0.init = vector("numeric", s)
  beta1.init = matrix(rep(0, s * p), s, p)
  for (sp in 1:s) {
    # tmp.m = stats::glm.fit(x = cbind(rep(1, nrow(X)), X), y = Y[,sp], family = family)
    tmp.m = glm2::glm.fit2(x = cbind(rep(1, nrow(X)), X), y = Y[,sp], family = family)
    beta0.init[sp] = tmp.m$coefficients[1]
    beta1.init[sp, ] = tmp.m$coefficients[2:(p + 1)]
  }
  clust = stats::kmeans(beta1.init, g)
  # update the starting parameters for the slopes and intercepts
  start.pars$beta0 = beta0.init
  start.pars$B = clust$centers
  if (!mixing.props.equal) {
    start.pars$pi = prop.table(table(clust$cluster))
  }
  attr(start.pars, "sp clust") = clust$clust

  if (update.intercepts) {
    for (sp in 1:s) {
      tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
      # tmp.m = stats::glm.fit(y = Y[,sp], x = matrix(1, nrow = nrow(X), ncol = 1), family = family, offset = tmp.offset)
      tmp.m = glm2::glm.fit2(y = Y[,sp], x = matrix(1, nrow = nrow(X), ncol = 1), family = family, offset = tmp.offset)
      start.pars$beta0[sp] = tmp.m$coefficients
    }
  }

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
#' @param fa.method Character string: indicating which FA software/approach to use.
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
init.fa.pars <- function(Y, X, g = 3, family = poisson(), d = 2, fa.method = c("svd", "prcomp", "qr", "fa", "glmmTMB", "gllvm"), update.intercepts = TRUE, constrain = TRUE, use.internal.qresid = FALSE, mixing.props.equal = FALSE) {

  fa.method = match.arg(fa.method)
  n <- nrow(Y); s <- ncol(Y); p <- ncol(X)

  # initialise usual SAM parameters
  start.pars = csam::init.sam.pars(Y, X, g = g, family = family, update.intercepts = update.intercepts, mixing.props.equal = mixing.props.equal)

  if (use.internal.qresid) {
    res <- csam:::qresiduals_csam_internal(start.pars, Y, X, family = family)
  } else {
    # fit an initial species-specific models to obtain warm starts as in Hui et al. 2013

    # with hard class labelling:
    res = matrix(rep(0, n * s), n, s)
    if (family$family == "binomial") {
      for (sp in 1:s) {
        if (update.intercepts) {
          tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
          tmp.m = stats::glm(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
          start.pars$beta0[sp] = tmp.m$coefficients
        } else {
          tmp.offset = start.pars$beta0[sp] + X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
          tmp.m = stats::glm(y ~ 0, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
        }
        res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m))
      }
      # remove boundaries that result in infinite residuals
      res[res == 1] = 1 -.Machine$double.eps
      res[res == 0] = .Machine$double.eps
      res = qnorm(res)
    } else {
      for (sp in 1:s) {
        tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
        tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
        if (update.intercepts) {
          start.pars$beta0[sp] = tmp.m$fit$par
        }
        res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
      }
    }
    # for (sp in 1:s) {
    #   if (family$family == "binomial") {
    #     if (update.intercepts) {
    #       tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
    #       tmp.m = stats::glm(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
    #       start.pars$beta0[sp] = tmp.m$coefficients
    #     } else {
    #       tmp.offset = start.pars$beta0[sp] + X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
    #       tmp.m = stats::glm(y ~ 0, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
    #     }
    #     res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m))
    #     # remove boundaries that result in infinite residuals
    #     res[res == 1] = 1 -.Machine$double.eps
    #     res[res == 0] = .Machine$double.eps
    #     res = qnorm(res)
    #   } else {
    #     tmp.offset = X %*% start.pars$B[attr(start.pars,"sp clust")[sp], ]
    #     tmp.m = glmmTMB::glmmTMB(y ~ 1, data = data.frame(y = Y[,sp]), family = family, offset = tmp.offset)
    #     if (update.intercepts) {
    #       start.pars$beta0[sp] = tmp.m$fit$par
    #     }
    #     res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth")
    #   }
    #   # res[ , sp] = glmmTMB:::residuals.glmmTMB(tmp.m, type = "dunn-smyth") # CON: uses an unexported function, throws error for binomial, requires glmmTMB() rather than glm() above
    #   # res[ , sp] = DHARMa:::residuals.DHARMa(DHARMa::simulateResiduals(tmp.m), quantileFunction = qnorm) # CON: uses a non-exported function
    #   # res[ , sp] = DHARMa::simulateResiduals(tmp.m)$scaledResiduals
    #   # res[ , sp] = statmod::qresiduals(tmp.m) # CON: leads to Inf/-Inf values
    # }

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
  }

  if (fa.method == "fa") {
    # perform factor analysis on the residuals
    tmp.cor = stats::cor(res)
    tmp.fa = psych::fa(tmp.cor, nfactors = d)
    tmp.fa = factanal(tmp.cor, factors = d)

    # update the starting parameters FA terms
    start.pars$U = unname(tmp.fa$scores)
    start.pars$Lambda = unname(suppressWarnings(as.data.frame(tmp.fa$loadings)[1:s, ]))

  } else if (fa.method == "glmmTMB") {

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

    # update the starting parameters FA terms
    start.pars$U = sapply(scores, cbind)
    start.pars$Lambda = tmp.fa$obj$env$report(tmp.fa$fit$parfull)$fact_load[[which(tmp.fa$modelInfo$grpVar == "site" & sapply(tmp.fa$modelInfo$reStruc$condReStruc, function(x){x$blockCode}) == 9)]]

  } else if (fa.method == "gllvm") {
    # perform factor analysis on the residuals
    tmp.fa = gllvm::gllvm(y = res, family = gaussian(), num.lv = d)

    # update the starting parameters FA terms
    start.pars$U = unname(tmp.fa$lvs)
    start.pars$Lambda = unname(coef(tmp.fa)$Species.scores)
  } else if (fa.method == "svd") {

    M <- svd(res, d, d)

    start.pars$Lambda <- M$v %*% diag(sqrt(M$d[1:d]))
    start.pars$U <- M$u %*% diag(1 / sqrt(M$d[1:d]))


  } else if (fa.method == "qr") {

    Q <- qr.Q(qr(res), complete = FALSE)
    start.pars$U <- Q[, 1:d]

    start.pars$Lambda <- crossprod(res, start.pars$U)

  } else if (fa.method == "prcomp") {

    pc <- prcomp(res,
      center = FALSE, # residuals already mean ~ 0
      scale. = FALSE,   # DO NOT rescale N(0,1) residuals
      rank. = d         # truncated PCA (efficient)
    )

    D_half <- diag(sqrt(pc$sdev[1:d]))

    start.pars$U <- pc$x %*% solve(D_half)
    start.pars$Lambda <- pc$rotation %*% D_half

  }

  # apply constraints to the factor analytic terms if desired
  if (constrain) {
    if (all(start.pars$U == 0)) {
      warning("no constraint applied to factor terms: all scores are zero")
    } else {
      try(assign("new.fa", correct.uv(start.pars$U, start.pars$Lambda)))
      if (exists("new.fa")) {
        start.pars$U <- new.fa$u
        start.pars$Lambda <- new.fa$v
      } else {
        warning("post-hoc constraint could not be applied")
      }
    }
  }

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

      if (family$family %in% c("binomial", "poisson")) {
        dev = NULL
      } else {
        dev <- family$dev.resids(y = yj, mu = mu, wt = rep(1, n))
      }
      ll_gj <- -0.5 * family$aic(y = yj, n = rep(1, n), mu = mu, wt = rep(1, n), dev = dev)
      comp[k] <- sum(ll_gj) + log(pi[k])
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
fitted_csam_internal <- function(par.list, X, tau = NULL, family,
                           type = c("response", "link"),
                           archetype_specific = FALSE) {

  type <- match.arg(type)

  if (!archetype_specific & is.null(tau)) {
    stop("tau must be supplied when computing posterior predicted means (i.e., when archetype_specific = FALSE)")
  }

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

#' Calculate the fitted values of a CSAM/SAM, \eqn{\mu_{ijk}}, and compute the posterior predictive mean (internal)
#'
#' @keywords internal
qresiduals_csam_internal <- function(par.list, Y, X, tau = NULL, family,
                                 type = c("dunn-smyth", "pearson"),
                                 archetype_specific = FALSE, seed = NULL) {

  type <- match.arg(type)

  # get the fitted values
  mu <- fitted_csam_internal(par.list = par.list, X = X, family = family, type = "response", archetype_specific = TRUE)

  # initialise the archetype specific cumulative distribution storage
  piG_ <- array(NA, dim = dim(mu))
  piG <- array(NA, dim = dim(mu))
  if (family$family == "binomial") {
    for (k in 1:dim(mu)[3]) {
      n <- matrix(1, nrow = nrow(Y), ncol = ncol(Y)) # may want to update if allowing different specifications for trials in a Bernoulli
      y <- n * Y
      if (!is.null(tau)) {

        # calculate the archetype CDFs and multiply by appropriate weights (tau[j, k], for each species/column j)
        piG_[ , , k] <- sweep(pbinom(y - 1, n, mu[ , , k]), 2, tau[ , k], "*")
        piG[ , , k] <- sweep(pbinom(y, n, mu[ , , k]), 2, tau[ , k], "*")

      } else {

        # calculate the archetype CDFs and multiply by appropriate weights (pi[k], for each archetype k)
        piG_[ , , k] <- pbinom(y - 1, n, mu[ , , k]) * par.list$pi[k]
        piG[ , , k] <- pbinom(y, n, mu[ , , k]) * par.list$pi[k]

      }
    }
  } else if (family$family == "poisson") {
    for (k in 1:dim(mu)[3]) {
      if (!is.null(tau)) {

        # calculate the archetype CDFs and multiply by appropriate weights (tau[j, k], for each species/column j)
        piG_[ , , k] <- sweep(ppois(Y - 1, n, mu[ , , k]), 2, tau[ , k], "*")
        piG[ , , k] <- sweep(ppois(Y, n, mu[ , , k]), 2, tau[ , k], "*")

      } else {

        # calculate the archetype CDFs and multiply by appropriate weights (pi[k], for each archetype k)
        piG_[ , , k] <- ppois(Y - 1, n, mu[ , , k]) * par.list$pi[k]
        piG[ , , k] <- ppois(Y, n, mu[ , , k]) * par.list$pi[k]

      }
    }
  } else {
    stop("family is not yet implemented")
  }

  # compute the sum of CDFs over the archetypes
  cf_ <- apply(piG_, c(1, 2), sum)
  cf <- apply(piG, c(1, 2), sum)
  # randomised standard normal
  u = runif(n = length(Y), min = cf_, max = cf)
  r = matrix(qnorm(u), nrow = nrow(Y), ncol = ncol(Y))

  return(r)
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

  loglik_mat <- matrix(0, s, g)

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

  ## ---- Add log pi and normalise ----
  logpost <- sweep(loglik_mat, 2, log(pi), "+")
  row_max <- apply(logpost, 1, max)
  tau_mat <- exp(logpost - row_max)
  tau <- tau_mat / rowSums(tau_mat)

  tau
}

#' standard error estimator for csam fits
#'
#' Computes standard errors for all model parameters in a fitted `"csam"` object using TMB's sdreport function.
#'
#' The fitted `csam` object **must** contain `Y`, `X`, and the penalty values
#' `psi1` and `psi2` (these are inferred from the object).
#'
#' @param object A fitted object of class `"csam"` (must contain `Y`, `X`,
#'   `psi1`, and `psi2`)..
#' @param U.random Logical; indicating whether to treat site effects as fixed or random when using \code{method = "tmb"}.
#' @param which.pars Character; of parameter group name(s) to be returned. Must be a subset of the canonical parameters of a CSAM TMB fit.
#'
#' @return A two column data frame of estimates and corresponding standard errors
#'
#' @examples
#' \dontrun{
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' fit$Y <- Y; fit$X <- X
#' # csam now stores psi1 and psi2 internally
#' V_louis <- vcov(fit, method = "louis")
#' }
#'
#' @export
#'
#' @importFrom TMB sdreport MakeADFun
se_csam <- function(object, U.random = TRUE, which.pars = NULL) {
  n <- nrow(object$Y)
  s <- ncol(object$Y)
  g <- nrow(object$B)

  # determine if the model is a vanilla SAM
  vanilla.sam <- is.null(object$Lambda)

  # set up required data list for the AD function
  data.list <- list(
    Y = object$Y,
    X = object$X,
    psi1 = if (!vanilla.sam) { object$psi1 } else { 0 },
    psi2 = if (!vanilla.sam) { object$psi2 } else { 0 },
    family = switch(object$family$family,
                    poisson = 0,
                    binomial = 1,
                    gaussian = 2),  # 0=Poisson, 1=Binomial, 2=Gaussian,...
    lik_type = as.numeric(U.random)
  )
  # set up required parameter list for the AD function
  start.pars.tmb = object$par.list
  start.pars.tmb$logit_pi = stats::binomial()$linkfun(object$par.list$pi[-g])
  start.pars.tmb$pi = NULL
  start.pars.tmb$log_phi = log(object$par.list$phi)
  start.pars.tmb$phi = NULL

  # create the mapping of parameters according to family
  mapping = list()
  if (object$family$family %in% c("binomial", "poisson")) {
    mapping$log_phi = factor(rep(NA, s))
  }
  # adjust inputs according to if the model is a vanilla SAM
  if (vanilla.sam) {
    start.pars.tmb$Lambda = matrix(rep(0, s), s, 1)
    start.pars.tmb$U <- matrix(rep(0, n), n, 1)
    data.list$lik_type = 0
    mapping$U = factor(rep(NA, n * 1))
    mapping$Lambda = factor(rep(NA, s * 1))
  }

  # set up the objective function
  if (U.random) {
    if (vanilla.sam) {
      warning("Marginalised likelihood for a vanilla SAM will be slightly off")
    }
    obj <- TMB::MakeADFun(
      data = data.list,
      parameters = start.pars.tmb,
      DLL = "csam", random = "U", map = mapping,
      silent = TRUE
    )
  } else {
    obj <- TMB::MakeADFun(
      data = data.list,
      parameters = start.pars.tmb,
      DLL = "csam", map = mapping,
      silent = TRUE
    )
  }

  out <- summary(TMB::sdreport(obj))

  if (!is.null(which.pars)) {
    if (!all(which.pars %in% c("beta0", "B", "logit_pi", "Lambda", "U", "log_phi"))) {
      stop(paste0(which.pars[!which.pars %in% c("beta0", "B", "logit_pi", "Lambda", "U", "log_phi")], " not found in canonical TMB parameters"))
    }
    out <- out[rownames(out) %in% which.pars, ]
  }

  return(out)
}

#' standard error using sandwich estimator for csam
#'
#' Computes standard errors for all model parameters in a fitted `"csam"` object using TMB's sdreport function.
#'
#' The fitted `csam` object **must** contain `Y`, `X`, and the penalty values
#' `psi1` and `psi2` (these are inferred from the object).
#'
#' @param object A fitted object of class `"csam"` (must contain `Y`, `X`,
#'   `psi1`, and `psi2`)..
#' @param U.random Logical; indicating whether to treat site effects as fixed or random when using \code{method = "tmb"}.
#' @param which.pars Character; of parameter group name(s) to be returned. Must be a subset of the canonical parameters of a CSAM TMB fit.
#'
#' @return A two column data frame of estimates and corresponding standard errors
#'
#' @examples
#' \dontrun{
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' fit$Y <- Y; fit$X <- X
#' # csam now stores psi1 and psi2 internally
#' V_louis <- vcov(fit, method = "louis")
#' }
#'
#' @export
#'
#' @importFrom TMB MakeADFun
#' @importFrom Matrix tcrossprod
#' @importFrom methods as
#' @importFrom MASS ginv
se_csam_sw <- function(object, which.pars = NULL) {
  n <- nrow(object$Y)
  s <- ncol(object$Y)
  g <- nrow(object$B)

  # determine if the model is a vanilla SAM
  vanilla.sam <- is.null(object$Lambda)

  # set up required data list for the AD function
  data.list_bread <- list(
    Y = object$Y,
    X = object$X,
    psi1 = if (!is.null(object$Lambda)) { object$psi1 } else { 0 },
    psi2 = if (!is.null(object$Lambda)) { object$psi2 } else { 0 },
    family = switch(object$family$family,
                    poisson = 0,
                    binomial = 1,
                    gaussian = 2),  # 0=Poisson, 1=Binomial, 2=Gaussian,...
    lik_type = 0
  )
  # set up required parameter list for the AD function
  par.list_bread = object$par.list
  par.list_bread$logit_pi = stats::binomial()$linkfun(object$par.list$pi[-g])
  par.list_bread$pi = NULL
  par.list_bread$log_phi = log(object$par.list$phi)
  par.list_bread$phi = NULL

  # create the mapping of parameters according to family
  mapping_bread = list()
  if (object$family$family %in% c("binomial", "poisson")) {
    mapping_bread$log_phi = factor(rep(NA, s))
  }
  # adjust inputs according to if the model is a vanilla SAM
  if (is.null(object$Lambda)) {
    par.list_bread$Lambda = matrix(rep(0, s), s, 1)
    par.list_bread$U <- matrix(rep(0, n), n, 1)
    mapping_bread$U = factor(rep(NA, n * 1))
    mapping_bread$Lambda = factor(rep(NA, s * 1))
  }

  obj_bread <- TMB::MakeADFun(
    data = data.list_bread,
    parameters = par.list_bread,
    DLL = "csam", map = mapping_bread,
    silent = TRUE
  )

  H = obj_bread$he()

  # Try and use a standard solve for inverting Hessian
  solve.type = 0
  try(assign("bread", solve(H)))
  # If that doesn't work try a robust approach
  # if (!exists("bread")) {
  #   # H <- (H + t(H)) / 2
  #   e <- eigen(H, symmetric = TRUE)
  #
  #   # Condition-based cutoff
  #   cutoff <- .Machine$double.eps * max(abs(e$values))
  #
  #   d_inv <- ifelse(abs(e$values) > cutoff, 1 / e$values, 0)
  #
  #   bread <- e$vectors %*% diag(d_inv) %*% t(e$vectors)
  #   solve.type = 1
  # }
  if (!exists("bread")) {
    bread <- MASS::ginv(H, tol = .Machine$double.eps)
    solve.type = 1
  }

  # get the indices of the parameters as appearing in the Hessian
  idx <- 1:length(obj_bread$par)
  names(idx) <- names(obj_bread$par)
  idx.list <- lapply(split(idx, names(idx)), unname)
  # adjust the parameters that occur in matrices
  idx.list$B <- matrix(idx.list$B, nrow = g)
  idx.list$Lambda <- matrix(idx.list$Lambda, nrow = s, ncol = d)
  idx.list$U <- matrix(idx.list$U, nrow = n, ncol = d)

  # work out the meat of the sandwich estimator

  # need to do this for each species
  S <- matrix(NA, nrow = nrow(bread), ncol = s)
  for (j in 1:s) {

    dat.list_meat <- list(
      Y = as.matrix(object$Y[ , j]),
      X =  as.matrix(object$X),
      psi1 = if (!is.null(object$Lambda)) { object$psi1 } else { 0 },
      psi2 = if (!is.null(object$Lambda)) { object$psi2 } else { 0 },
      family = switch(object$family$family,
                      poisson = 0,
                      binomial = 1,
                      gaussian = 2),
      lik_type = 0
    )

    par.list_meat = list(
      beta0    = object$beta0[j],
      B        = object$B,
      logit_pi = binomial()$linkfun(object$pi[1:(g - 1)]),
      U        = object$U,
      Lambda   = matrix(object$Lambda[j, ], nrow = 1),
      log_phi  = object$phi[j]  # only used for Gaussian/Gamma/NB
    )

    # create the mapping of parameters according to family
    mapping_meat = list()
    if (object$family$family %in% c("binomial", "poisson")) {
      mapping_meat$log_phi = factor(rep(NA, 1))
    }

    obj_meat <- TMB::MakeADFun(
      data = dat.list_meat,
      parameters = par.list_meat,
      DLL = "csam", map = mapping_meat,
      silent = TRUE
    )

    S_j <- as.vector(-obj_meat$gr())

    tmp.idx <- c(idx.list$beta0[j], as.vector(idx.list$B), idx.list$logit_pi, as.vector(idx.list$U), as.vector(idx.list$Lambda[j, ]), idx.list$logphi)
    S[tmp.idx, j] <- S_j
    rm(tmp.idx, dat.list_meat, par.list_meat, mapping_meat, obj_meat)
  }

  # put zeroes in for the NA values
  S[is.na(S)] <- 0
  # make sparse before computing meat
  S <- methods::as(S, "sparseMatrix")

  # compute the meat
  meat = Matrix::tcrossprod(S)

  # compute the covariance matrix
  cov.mat <- bread %*% meat %*% bread
  se <- sqrt(diag(as.matrix(cov.mat)))

  # also compute the original standard errors for reference
  se_og <-  sqrt(diag(as.matrix(bread)))

  if (!is.null(which.pars)) {
    if (!all(which.pars %in% c("beta0", "B", "logit_pi", "Lambda", "U", "log_phi"))) {
      stop(paste0(which.pars[!which.pars %in% c("beta0", "B", "logit_pi", "Lambda", "U", "log_phi")], " not found in canonical TMB parameters"))
    }
    se <- se[names(obj_bread$par) == which.pars]
    se_og <- se_og[names(obj_bread$par) == which.pars]
  }
  attr(se, "solve.type") = solve.type
  attr(se, "og") = se_og

  return(se)
}

#' Impose constraints on the factor scores and loadings
#'
#' @description
#' Imposes constraints on the factor scores and loadings as in GMF.
#'
#' @param U A numeric matrix of factor scores of dimension \eqn{n \times d},
#'   with rows corresponding to sites and columns to the dimension of the factor analytic term.
#' @param V A numeric matrix of factor loadings of dimension \eqn{s \times d},
#'   with rows corresponding to species and columns to the dimension of the factor analytic term.
#'
#' @return A list of the constrained factor analytic terms with components:
#' \describe{
#'   \item{u}{A numeric matrix of factor scores of dimension \eqn{n \times d}}
#'   \item{v}{A numeric matrix of factor loadings of dimension \eqn{s \times d}}
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
#' @importFrom whitening whiteningMatrix
correct.uv <- function (U, V) {
  S = cov(U)
  if (ncol(U) == 1) {
    return(list(u = U/sqrt(c(S)), v = V * sqrt(c(S))))
  }
  W = whitening::whiteningMatrix(S)
  U = U %*% W
  V = V %*% t(solve(W))
  V.qr = qr(t(V))
  U = U %*% qr.Q(V.qr)
  V = t(qr.R(V.qr))
  d = diag(V)
  V = t(sign(d) * t(V))
  U = t(sign(d) * t(U))
  list(u = U, v = V)
}
