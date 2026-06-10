##### CSAM EM algorithm components #####

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
  XB <- tcrossprod(X, B)          # n × g
  UL <- tcrossprod(U, Lambda)     # n × s

  loglik_mat <- matrix(0, s, g)

  ## ---- E-step loops (species-level, parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    eta_base <- beta0[j] + UL[, j]

    for (k in seq_len(g)) {

      eta <- eta_base + XB[, k]
      mu  <- family$linkinv(eta)

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

  ## ---- Posterior normalisation (softmax) ----
  logpost <- sweep(loglik_mat, 2, log(pi), "+")

  ## stabilise row-wise
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
                            family, maxit = 1, use.starts = FALSE) {

  beta0 = par.list$beta0; B = par.list$B; pi = par.list$pi; U = par.list$U; Lambda = par.list$Lambda; phi = par.list$phi

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X)
  B_new <- B

  ## ---- Shared linear predictor (fixed in this CM-step) ----
  UL <- tcrossprod(U, Lambda)

  ## ---- Pre-build stacked design matrix ONCE ----
  # Stack X by species
  X_stack <- X[rep(seq_len(n), times = s), ]

  # compute the offsets which are constant for each archetype
  offset_stack <- rep(beta0, each = n) + as.vector(UL)

  # put the response data in long format
  y_stack <- as.vector(Y)

  for (k in seq_len(g)) {

    weights_stack <- rep(tau[ , k], each = n)

    # fit <- stats::glm.fit(
    #   x = X_stack,
    #   y = y_stack,
    #   weights = weights_stack,
    #   start = if (use.starts) {B_new[k, ]} else {NULL},
    #   offset = offset_stack,
    #   family = family,
    #   control = list(maxit = maxit)
    # )
    fit <- glm2::glm.fit2(
      x = X_stack,
      y = y_stack,
      weights = weights_stack,
      start = if (use.starts) {B_new[k, ]} else {NULL},
      offset = offset_stack,
      family = family,
      control = list(maxit = maxit)
    )

    if (exists("fit")) {
      if (is.null(fit$coefficients) || any(is.na(fit$coefficients)) || !all(is.finite(fit$coefficients))) {
        warning(sprintf("glm.fit failed for archetype %d; leaving B[k,] unchanged", k))
      } else {
        B_new[k, ] <- fit$coefficients
      }
    } else {
      warning(sprintf("glm.fit failed for archetype %d; leaving B[k,] unchanged", k))
    }
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
    family, psi2 = 0, maxit = 1, backend = c("C++", "R"), use.starts = FALSE, use.glm.fit.when.unpenalised = FALSE) {
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
  XB <- tcrossprod(X, B)          # n × g
  XU <- cbind(1, U)         # n × (d+1)

  beta0_new  <- beta0
  Lambda_new <- Lambda
  phi_new    <- phi

  # compute the design matrix and offsets (these are constant across species given archetypes)
  X_stack <- XU[rep(seq_len(nrow(XU)), g), ]
  offset_stack <- as.vector(XB)

  ## ---- Species loop (parallelisable) ----
  for (j in seq_len(s)) {

    yj <- Y[, j]
    y_stack <- rep(Y[ , j], g)
    weights_stack <- rep(tau[j, ], each = n)

    if (psi2 == 0 & use.glm.fit.when.unpenalised) {
      # fit <- stats::glm.fit(
      #   x = X_stack,
      #   y = y_stack,
      #   weights = weights_stack,
      #   start = if (use.starts) {c(beta0_new[j], Lambda_new[j, ])} else {NULL},
      #   offset = offset_stack,
      #   family = family,
      #   control = list(maxit = maxit)
      # )
      fit <- glm2::glm.fit2(
        x = X_stack,
        y = y_stack,
        weights = weights_stack,
        start = if (use.starts) {c(beta0_new[j], Lambda_new[j, ])} else {NULL},
        offset = offset_stack,
        family = family,
        control = list(maxit = maxit)
      )
    } else {
      fit <- penalised_glm_fit(
        X       = X_stack,
        y       = y_stack,
        weights = weights_stack,
        offset  = offset_stack,
        penalty_weights = c(0, rep(1, d)),
        lambda  = psi2,
        family  = family,
        start   = if (use.starts) {c(beta0_new[j], Lambda_new[j, ])} else {NULL},
        maxit   = maxit,
        backend = backend
      )
    }

    if (exists("fit")) {
      if (is.null(fit$coefficients) || any(is.na(fit$coefficients)) || !all(is.finite(fit$coefficients))) {
        warning(sprintf("(P-)IRLS failed for species %d; leaving parameters unchanged", j))
      } else {
        beta0_new[j]    <- fit$coefficients[1]
        Lambda_new[j, ] <- fit$coefficients[-1]
        if (family$family == "gaussian") {
          phi_new[j] <- sum(weights_stack * (y_stack - fit$mu)^2) /
            sum(weights_stack)
        }
      }
    } else {
      warning(sprintf("(P-)IRLS failed for species %d; leaving parameters unchanged", j))
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
                              family, psi1 = 0, maxit = 1, backend = c("C++", "R"), use.starts = FALSE, use.glm.fit.when.unpenalised = FALSE) {
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
  XB <- tcrossprod(X, B)   # n × g

  U_new <- U

  # compute the design matrix and weights (these are constant over sites given archetypes)
  X_stack <- Lambda[rep(seq_len(nrow(Lambda)), g), ]
  weights_stack <- as.vector(tau)

  penalty_weights <- rep(1, d)

  for (i in seq_len(n)) {

    y_stack <- rep(Y[i, ], g)
    offset_stack <- rep(beta0, g) + rep(XB[i, ], each = s)

    if (psi1 == 0 & use.glm.fit.when.unpenalised) {
      # fit <- stats::glm.fit(
      #   x = X_stack,
      #   y = y_stack,
      #   weights = weights_stack,
      #   start = if (use.starts) {U_new[i, ]} else {NULL},
      #   offset = offset_stack,
      #   family = family,
      #   control = list(maxit = maxit)
      # )
      fit <- glm2::glm.fit2(
        x = X_stack,
        y = y_stack,
        weights = weights_stack,
        start = if (use.starts) {U_new[i, ]} else {NULL},
        offset = offset_stack,
        family = family,
        control = list(maxit = maxit)
      )
    } else {
      fit <- penalised_glm_fit(
        X = X_stack,
        y = y_stack,
        weights = weights_stack,
        offset = offset_stack,
        penalty_weights = penalty_weights,
        lambda = psi1,
        family = family,
        start = if (use.starts) {U_new[i, ]} else {NULL},
        maxit = maxit,
        backend = backend
      )
    }

    if (exists("fit")) {
      if (is.null(fit$coefficients) || any(is.na(fit$coefficients)) || !all(is.finite(fit$coefficients))) {
        warning(sprintf("(P-)IRLS failed for site %d; leaving U[i, ] unchanged", i))
      } else {
        U_new[i, ] <- fit$coefficients
      }
    } else {
      warning(sprintf("(P-)IRLS failed for site %d; leaving U[i, ] unchanged", i))
    }

  }

  U_new
}
