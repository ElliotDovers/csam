#' Plot diagnostics for CSAM fits
#'
#' @description
#' S3 method for the base `graphics::plot` generic that produces ECM traces
#' and diagnostic plots for objects of class `"csam"`. Available `param`
#' options are:
#' \itemize{
#'   \item `"pll"`: penalised log-likelihood trace,
#'   \item `"B"`: archetype coefficient traces (one panel per predictor),
#'   \item `"beta0"`: species intercept traces,
#'   \item `"pi"`: mixing proportion traces,
#'   \item `"U"`: site-score traces (one panel per factor),
#'   \item `"Lambda"`: species-loading traces (one panel per factor),
#'   \item `"phi"`: species dispersion traces (Gaussian family only),
#'   \item `"cor"`: correlation matrix induced by the estimated loadings
#'     (requires `corrplot` for a prettier display).
#' }
#'
#' @param x A fitted object of class `"csam"`.
#' @param param Character; which quantity to plot. See Details.
#' @param truth Optional list of true parameter values (named like the
#'   components of a `"csam"` object) to overlay as horizontal reference
#'   lines or to show the true correlation matrix in the upper triangle.
#' @param ylim Optional numeric vector of length 2 giving y-axis limits.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This function expects the fitted object to contain a `trace` element
#' produced when the model was fit with `trace = TRUE`. The trace should
#' contain lists of parameter values across ECM iterations (see the
#' package documentation for the structure).
#'
#' @examples
#' \dontrun{
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' plot(fit, param = "pll")
#' plot(fit, param = "cor")
#' }
#'
#' @exportS3Method graphics::plot csam
plot.csam <- function(x,
                      param = c("pll", "B", "beta0", "pi", "U", "Lambda", "phi", "cor"),
                      truth = NULL,
                      ylim = NULL,
                      ...) {

  if (is.null(x$trace)) {
    stop("No trace stored in fit object. Re-run csam() with trace = TRUE.")
  }

  param <- match.arg(param)
  tr <- x$trace

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Helper: distinct colours
  colfun <- function(n) grDevices::rainbow(n, s = 0.7, v = 0.9)

  # Helper: compute y-limits including truth if provided
  compute_ylim <- function(mat, truth_vals = NULL, user_ylim = NULL) {
    if (!is.null(user_ylim)) return(user_ylim)
    if (!is.null(truth_vals)) return(range(c(mat, truth_vals)))
    range(mat)
  }

  # -----------------------------
  # 1. Penalised log-likelihood
  # -----------------------------
  if (param == "pll") {
    ylims <- compute_ylim(tr$pll, truth$pll %||% NULL, ylim)
    plot(tr$pll, type = "l", lwd = 2, col = "steelblue",
         xlab = "Iteration", ylab = "Penalised log-likelihood",
         main = "ECM Trace: Penalised Log-Likelihood",
         ylim = ylims)
    if (!is.null(truth) && !is.null(truth$pll)) {
      abline(h = truth$pll, lty = 2, col = "black")
    }
    return(invisible())
  }

  # -----------------------------
  # 2. B (g x p) — one panel per predictor
  # -----------------------------
  if (param == "B") {
    B_list <- tr$B
    g <- nrow(x$B)
    p <- ncol(x$B)

    B_traces <- vector("list", p)
    for (j in 1:p) {
      B_traces[[j]] <- t(sapply(B_list, function(B) B[, j]))
    }

    cols <- colfun(g)

    op <- par(mfrow = c(p, 1), mar = c(4, 4, 3, 1))
    on.exit(par(op), add = TRUE)

    for (j in 1:p) {
      mat <- B_traces[[j]]
      truth_vals <- if (!is.null(truth) && !is.null(truth$B)) truth$B[, j] else NULL
      ylims <- compute_ylim(mat, truth_vals, ylim)

      plot(mat[, 1], type = "l", col = cols[1], lwd = 2,
           xlab = "Iteration", ylab = paste0("B[,", j, "]"),
           main = paste0("ECM Trace: B predictor ", j),
           ylim = ylims)

      for (k in 2:g) lines(mat[, k], col = cols[k], lwd = 2)

      if (!is.null(truth_vals)) {
        for (k in 1:g) abline(h = truth_vals[k], col = cols[k], lty = 2)
      }

      legend("topright",
             legend = paste0("B[", 1:g, ",", j, "]"),
             col = cols, lwd = 2, cex = 0.7)
    }

    return(invisible())
  }

  # -----------------------------
  # 3. beta0
  # -----------------------------
  if (param == "beta0") {
    beta_list <- tr$beta0
    s <- length(x$beta0)

    mat <- t(sapply(beta_list, identity))
    cols <- colfun(s)

    truth_vals <- if (!is.null(truth)) truth$beta0 else NULL
    ylims <- compute_ylim(mat, truth_vals, ylim)

    plot(mat[, 1], type = "l", col = cols[1], lwd = 2,
         xlab = "Iteration", ylab = "beta0",
         main = "ECM Trace: beta0",
         ylim = ylims)

    for (j in 2:s) lines(mat[, j], col = cols[j], lwd = 2)

    if (!is.null(truth_vals)) {
      for (j in 1:s) abline(h = truth_vals[j], col = cols[j], lty = 2)
    }

    legend("topright", legend = paste0("beta0[", 1:s, "]"),
           col = cols, lwd = 2, cex = 0.7)
    return(invisible())
  }

  # -----------------------------
  # 4. pi
  # -----------------------------
  if (param == "pi") {
    pi_list <- tr$pi
    g <- length(x$pi)

    mat <- t(sapply(pi_list, identity))
    cols <- colfun(g)

    truth_vals <- if (!is.null(truth)) truth$pi else NULL
    ylims <- compute_ylim(mat, truth_vals, ylim)

    plot(mat[, 1], type = "l", col = cols[1], lwd = 2,
         xlab = "Iteration", ylab = "pi",
         main = "ECM Trace: Mixing Proportions",
         ylim = ylims)

    for (k in 2:g) lines(mat[, k], col = cols[k], lwd = 2)

    if (!is.null(truth_vals)) {
      for (k in 1:g) abline(h = truth_vals[k], col = cols[k], lty = 2)
    }

    legend("topright", legend = paste0("pi[", 1:g, "]"),
           col = cols, lwd = 2, cex = 0.7)
    return(invisible())
  }

  # -----------------------------
  # 5. U (n x d)
  # -----------------------------
  if (param == "U") {
    U_list <- tr$U
    n <- nrow(x$U)
    d <- ncol(x$U)

    op <- par(mfrow = c(d, 1), mar = c(4, 4, 3, 1))
    on.exit(par(op), add = TRUE)

    for (h in 1:d) {
      mat_h <- t(sapply(U_list, function(U) U[, h]))
      cols <- colfun(n)

      truth_vals <- if (!is.null(truth)) truth$U[, h] else NULL
      ylims <- compute_ylim(mat_h, truth_vals, ylim)

      plot(mat_h[, 1], type = "l", col = cols[1], lwd = 1,
           xlab = "Iteration", ylab = paste0("U[,", h, "]"),
           main = paste0("ECM Trace: U factor ", h),
           ylim = ylims)

      for (i in 2:n) lines(mat_h[, i], col = cols[i], lwd = 1)

      if (!is.null(truth_vals)) {
        for (i in 1:n) abline(h = truth_vals[i], col = cols[i], lty = 2)
      }

      legend("topright",
             legend = paste0("U[", 1:n, ",", h, "]"),
             col = cols, lwd = 1, cex = 0.5)
    }

    return(invisible())
  }

  # -----------------------------
  # 6. Lambda (s x d)
  # -----------------------------
  if (param == "Lambda") {
    L_list <- tr$Lambda
    s <- nrow(x$Lambda)
    d <- ncol(x$Lambda)

    op <- par(mfrow = c(d, 1), mar = c(4, 4, 3, 1))
    on.exit(par(op), add = TRUE)

    for (h in 1:d) {
      mat_h <- t(sapply(L_list, function(L) L[, h]))
      cols <- colfun(s)

      truth_vals <- if (!is.null(truth)) truth$Lambda[, h] else NULL
      ylims <- compute_ylim(mat_h, truth_vals, ylim)

      plot(mat_h[, 1], type = "l", col = cols[1], lwd = 1,
           xlab = "Iteration", ylab = paste0("Lambda[,", h, "]"),
           main = paste0("ECM Trace: Lambda factor ", h),
           ylim = ylims)

      for (j in 2:s) lines(mat_h[, j], col = cols[j], lwd = 1)

      if (!is.null(truth_vals)) {
        for (j in 1:s) abline(h = truth_vals[j], col = cols[j], lty = 2)
      }

      legend("topright",
             legend = paste0("Lambda[", 1:s, ",", h, "]"),
             col = cols, lwd = 1, cex = 0.5)
    }

    return(invisible())
  }

  # -----------------------------
  # 7. phi
  # -----------------------------
  if (param == "phi") {
    phi_list <- tr$phi
    s <- length(x$phi)

    mat <- t(sapply(phi_list, identity))
    cols <- colfun(s)

    truth_vals <- if (!is.null(truth)) truth$phi else NULL
    ylims <- compute_ylim(mat, truth_vals, ylim)

    plot(mat[, 1], type = "l", col = cols[1], lwd = 2,
         xlab = "Iteration", ylab = "phi",
         main = "ECM Trace: phi parameters",
         ylim = ylims)

    for (j in 2:s) lines(mat[, j], col = cols[j], lwd = 2)

    if (!is.null(truth_vals)) {
      for (j in 1:s) abline(h = truth_vals[j], col = cols[j], lty = 2)
    }

    legend("topright", legend = paste0("phi[", 1:s, "]"),
           col = cols, lwd = 2, cex = 0.7)
    return(invisible())
  }

  # -----------------------------
  # 8. Correlation matrix induced by Lambda
  # -----------------------------
  if (param == "cor") {

    if (!requireNamespace("corrplot", quietly = TRUE)) {
      stop("Package 'corrplot' is required for correlation plotting.")
    }

    Lambda_est <- x$Lambda
    cor_est <- stats::cov2cor(Lambda_est %*% t(Lambda_est))

    truth_supplied <- !is.null(truth) && !is.null(truth$Lambda)

    if (truth_supplied) {
      Lambda_true <- truth$Lambda
      cor_true <- stats::cov2cor(Lambda_true %*% t(Lambda_true))
      cor_est[upper.tri(cor_est)] <- cor_true[upper.tri(cor_true)]
    }

    corrplot::corrplot(
      cor_est,
      method = "square",
      tl.cex = 1,
      tl.srt = 45,
      diag = TRUE,
      tl.col = "black",
      addCoef.col = "black",
      mar = c(0, 2.1, 2.1, 0)
    )

    if (truth_supplied) {
      S <- nrow(cor_est)
      abline(a = S + 1, b = -1, xpd = TRUE, lty = "dashed")

      mtext("Estimated", side = 2, line = 3, at = 3)
      mtext("Truth", side = 3, line = 2.5, at = 4.25)
    }

    return(invisible())
  }

  stop("Unknown parameter type.")
}

#' Variance-covariance estimators for csam fits (analytic, modular)
#'
#' S3 method for the base `stats::vcov` generic that computes variance-covariance
#' matrices for all model parameters in a fitted `"csam"` object using analytic
#' complete-data scores and expected Hessians. Supported `method` values are
#' `"naive"`, `"louis"`, `"oakes"`, and `"sandwich"`.
#'
#' The returned matrix is on an unconstrained parameter scale:
#' - **beta0**: raw (length s)
#' - **B**: stacked by archetype then predictor (length g * p)
#' - **pi**: log-odds for first (g-1) components vs component g
#' - **U**: stacked site scores (length n * d)
#' - **Lambda**: stacked loadings (length s * d)
#' - **phi**: log(phi) (length s) for Gaussian family
#'
#' The fitted `csam` object **must** contain `Y`, `X`, and the penalty values
#' `psi1` and `psi2` (these are inferred from the object and therefore not
#' required as arguments to `vcov`).
#'
#' @param object A fitted object of class `"csam"` (must contain `Y`, `X`,
#'   `psi1`, and `psi2`).
#' @param method Character; one of `"naive"`, `"louis"`, `"oakes"`, `"sandwich"`.
#' @param ... Currently unused.
#'
#' @return A symmetric variance-covariance matrix for the full parameter
#'   vector (on the unconstrained scale). Row/column names indicate the
#'   parameter ordering.
#'
#' @examples
#' \dontrun{
#' fit <- csam(Y, X, g = 2, d = 1, family = poisson(), trace = TRUE)
#' fit$Y <- Y; fit$X <- X
#' # csam now stores psi1 and psi2 internally
#' V_louis <- vcov(fit, method = "louis")
#' }
#'
#' @exportS3Method stats::vcov csam
vcov.csam <- function(object,
                       method = c("naive", "louis", "tmb", "sandwich"),
                       ...) {

  method <- match.arg(method)

  if (method == "naive") {

    S <- score_complete(object)        # (s*g) × pθ
    H <- hessian_complete(object)      # list of length s*g
    tau <- object$tau                  # s × g

    s <- ncol(object$Y)
    g <- nrow(object$B)

    p_theta <- ncol(S)
    name_vec <- colnames(S)

    tau_vec <- as.vector(t(tau))       # length s*g, row-major (j,k)
    Info <- Reduce(`+`, Map(function(Hjk, w) w * Hjk, H, tau_vec))
    covmat <- tryCatch(solve(Info), error = function(e) MASS::ginv(Info))
    colnames(covmat) <- rownames(covmat) <- name_vec

  } else if (method == "louis") {

    S <- score_complete(object)        # (s*g) × pθ
    H <- hessian_complete(object)      # list of length s*g
    tau <- object$tau                  # s × g

    s <- ncol(object$Y)
    g <- nrow(object$B)

    p_theta <- ncol(S)
    name_vec <- colnames(S)

    tau_vec <- as.vector(t(tau))       # length s*g, row-major (j,k)

    Info <- matrix(0, p_theta, p_theta)

    for (j in 1:s) {

      idx_j <- ((j-1)*g + 1):(j*g)

      Hbar_j <- Reduce(`+`, Map(function(Hjk, w) w * Hjk, H[idx_j], tau[j, ]))

      S_j <- S[idx_j, , drop = FALSE]  # g × p_theta

      sbar_j <- colSums(tau[j, ] * S_j)

      S2_j <- matrix(0, p_theta, p_theta)
      for (k in 1:g) {
        sjk <- S_j[k, ]
        S2_j <- S2_j + tau[j, k] * tcrossprod(sjk)
      }

      Info <- Info + Hbar_j - (S2_j - tcrossprod(sbar_j))
    }

    covmat <- tryCatch(solve(Info), error = function(e) MASS::ginv(Info))
    colnames(covmat) <- rownames(covmat) <- name_vec

  } else if (method == "tmb") {

    # set up required data list for the AD function
    data_list <- list(
      Y = object$Y,
      X = object$X,
      psi1 = object$psi1,
      psi2 = object$psi2,
      family = switch(object$family$family,
                      poisson = 0,
                      binomial = 1,
                      gaussian = 2),  # 0=Poisson, 1=Binomial, 2=Gaussian,...
      lik_type = 1
    )
    # set up required parameter list for the AD function
    start.pars.tmb = m[1:6]
    start.pars.tmb$theta_pi = log(m$pi)
    start.pars.tmb$pi = NULL
    start.pars.tmb$logphi = log(m$phi)
    start.pars.tmb$phi = NULL
    obj <- TMB::MakeADFun(
      data = data_list,
      parameters = start.pars.tmb,
      DLL = "csam", random = "U", map = list(logphi = factor(rep(NA, s)), theta_pi = factor(rep(NA, g))), # for now fixing these to cheat into non-singular matrix
      silent = TRUE
    )
    rep <- TMB::sdreport(obj, getJointPrecision = TRUE)
    covmat = solve(as.matrix(rep$jointPrecision))

  } else if (method == "sandwich") {

    S <- score_complete(object)        # (s*g) × pθ
    H <- hessian_complete(object)      # list of length s*g
    tau <- object$tau                  # s × g

    s <- ncol(object$Y)
    g <- nrow(object$B)

    p_theta <- ncol(S)
    name_vec <- colnames(S)

    tau_vec <- as.vector(t(tau))       # length s*g, row-major (j,k)

    Info <- matrix(0, p_theta, p_theta)
    Meat <- matrix(0, p_theta, p_theta)

    for (j in 1:s) {

      idx_j <- ((j-1)*g + 1):(j*g)

      Hbar_j <- Reduce(`+`, Map(function(Hjk, w) w * Hjk, H[idx_j], tau[j, ]))

      S_j <- S[idx_j, , drop = FALSE]  # g × p_theta

      sbar_j <- colSums(tau[j, ] * S_j)

      S2_j <- matrix(0, p_theta, p_theta)
      for (k in 1:g) {
        sjk <- S_j[k, ]
        S2_j <- S2_j + tau[j, k] * tcrossprod(sjk)
      }

      Info <- Info + Hbar_j - (S2_j - tcrossprod(sbar_j))

      Meat <- Meat + tcrossprod(sbar_j)
    }

    bread <- tryCatch(solve(Info), error = function(e) MASS::ginv(Info))
    covmat <- bread %*% Meat %*% bread
    colnames(covmat) <- rownames(covmat) <- name_vec

  } else {
    stop("Unknown method")
  }

  return(covmat)
}

#' Confidence intervals for csam fits (simple vcov-based, excluding pi and phi)
#'
#' @description
#' Compute Wald confidence intervals using standard errors from vcov(object).
#' This method excludes mixing proportions (pi) and dispersion (phi) from
#' interval calculation. Use `which` to select a block: "beta0", "B", "U", "Lambda".
#'
#' @exportS3Method stats::confint csam
confint.csam <- function(object,
                         method = c("naive", "louis", "tmb", "sandwich"),
                         level = 0.95,
                         which = NULL,
                         ...) {

  method <- match.arg(method)

  if (!inherits(object, "csam"))
    stop("object must be of class 'csam'")

  Y <- object$Y
  X <- object$X
  family <- object$family

  s <- ncol(Y)
  g <- nrow(object$B)
  p <- ncol(X)
  n <- nrow(object$U)
  d <- ncol(object$U)
  is_gaussian <- identical(family$family, "gaussian")

  # ------------------------------------------------------------
  # 1. Extract parameter vector in canonical order
  # ------------------------------------------------------------
  t_pi <- if (g > 1) log(object$pi[1:(g-1)] / object$pi[g]) else numeric(0)
  logphi <- if (is_gaussian) log(object$phi) else numeric(0)

  theta <- c(
    object$beta0,
    as.vector(t(object$B)),
    t_pi,
    as.vector(t(object$U)),
    as.vector(t(object$Lambda)),
    logphi
  )

  # ------------------------------------------------------------
  # 2. Parameter names (must match vcov.csam)
  # ------------------------------------------------------------
  param_names <- c(
    paste0("beta0[", seq_len(s), "]"),
    unlist(lapply(seq_len(g), function(k)
      paste0("B[", k, ",", seq_len(p), "]"))),
    if (g > 1) paste0("t_pi[", seq_len(g - 1), "]") else character(0),
    unlist(lapply(seq_len(n), function(i)
      paste0("U[", i, ",", seq_len(d), "]"))),
    unlist(lapply(seq_len(s), function(j)
      paste0("Lambda[", j, ",", seq_len(d), "]"))),
    if (is_gaussian) paste0("logphi[", seq_len(s), "]") else character(0)
  )

  names(theta) <- param_names

  # ------------------------------------------------------------
  # 3. Covariance matrix and standard errors
  # ------------------------------------------------------------
  V <- vcov.csam(object, method = method, ...)
  se <- sqrt(diag(V))

  # ------------------------------------------------------------
  # 4. Confidence intervals
  # ------------------------------------------------------------
  alpha <- 1 - level
  z <- qnorm(1 - alpha/2)

  lower <- theta - z * se
  upper <- theta + z * se

  # ------------------------------------------------------------
  # 5. Parameter selection (optional)
  # ------------------------------------------------------------
  if (!is.null(which)) {

    if (which %in% param_names) {
      idx <- which(param_names == which)
    } else {
      idx <- grepl(which, param_names, fixed = TRUE)
    }

    theta <- theta[idx]
    lower <- lower[idx]
    upper <- upper[idx]
    param_names <- param_names[idx]
  }

  # ------------------------------------------------------------
  # 6. Return tidy data frame
  # ------------------------------------------------------------
  ret <- data.frame(lower, upper)
  colnames(ret) <- c(
    paste0(100 * alpha/2, "%"),
    paste0(100 * (1 - alpha/2), "%")
  )
  rownames(ret) <- param_names

  ret
}
