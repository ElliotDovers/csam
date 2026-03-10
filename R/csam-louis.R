score_all <- function(object) {

  if (!inherits(object, "csam"))
    stop("object must be a fitted csam model")

  Y      <- object$Y
  X      <- object$X
  U      <- object$U
  beta0  <- object$beta0
  B      <- object$B
  Lambda <- object$Lambda
  phi    <- object$phi
  pi     <- object$pi
  family <- object$family

  n <- nrow(Y)
  s <- ncol(Y)
  p <- ncol(X)
  g <- nrow(B)
  d <- ncol(U)

  ## a(φ) for canonical families
  a_fun <- function(fam, phi_j) {
    switch(fam,
           "poisson"  = 1,
           "binomial" = 1,
           "Gamma"    = phi_j,
           stop("score_all: a(φ) not implemented for family ", fam)
    )
  }

  ## allocate arrays
  score_beta0  <- array(0, dim = c(s, g, 1))
  score_B      <- array(0, dim = c(s, g, p))
  score_Lambda <- array(0, dim = c(s, g, d))
  score_phi    <- array(0, dim = c(s, g, 1))
  score_pi     <- array(0, dim = c(s, g, 1))

  fam_name <- family$family

  for (j in 1:s) {
    yj       <- Y[, j]
    beta0_j  <- beta0[j]
    lambda_j <- Lambda[j, ]
    phi_j    <- phi[j]
    a_phi    <- a_fun(fam_name, phi_j)

    for (k in 1:g) {

      beta_k <- B[k, ]

      ## η_{ijk}
      eta <- as.vector(beta0_j + X %*% beta_k + U %*% lambda_j)

      ## μ_{ijk}
      mu  <- family$linkinv(eta)

      ## μ'(η)
      mu_eta <- family$mu.eta(eta)

      ## ∂ℓ/∂η = (y - μ)/(a(φ) * μ'(η))
      d_eta <- (yj - mu) / (a_phi * mu_eta)

      ## β0_j
      score_beta0[j, k, 1] <- sum(d_eta)

      ## β_k
      score_B[j, k, ] <- as.vector(t(X) %*% d_eta)

      ## λ_j
      score_Lambda[j, k, ] <- as.vector(t(U) %*% d_eta)

      ## φ_j (zero for Poisson/Binomial)
      if (fam_name %in% c("poisson", "binomial")) {
        score_phi[j, k, 1] <- 0
      } else {
        stop("score for φ_j not implemented for this family")
      }

      ## π_k
      score_pi[j, k, 1] <- 1 / pi[k]
    }
  }

  list(
    beta0   = score_beta0,
    B       = score_B,
    Lambda  = score_Lambda,
    phi     = score_phi,
    pi      = score_pi
  )
}

hessian_all <- function(object) {

  if (!inherits(object, "csam"))
    stop("object must be a fitted csam model")

  Y      <- object$Y
  X      <- object$X
  U      <- object$U
  beta0  <- object$beta0
  B      <- object$B
  Lambda <- object$Lambda
  phi    <- object$phi
  pi     <- object$pi
  family <- object$family

  n <- nrow(Y)
  s <- ncol(Y)
  p <- ncol(X)
  g <- nrow(B)
  d <- ncol(U)

  ## a(φ) for canonical families
  a_fun <- function(fam, phi_j) {
    switch(fam,
           "poisson"  = 1,
           "binomial" = 1,
           "Gamma"    = phi_j,
           stop("hessian_all: a(φ) not implemented for family ", fam)
    )
  }

  ## allocate arrays
  H_beta0_beta0   <- array(0, dim = c(s, g, 1, 1))
  H_B_B           <- array(0, dim = c(s, g, p, p))
  H_Lambda_Lambda <- array(0, dim = c(s, g, d, d))
  H_beta0_B       <- array(0, dim = c(s, g, 1, p))
  H_beta0_Lambda  <- array(0, dim = c(s, g, 1, d))
  H_B_Lambda      <- array(0, dim = c(s, g, p, d))
  H_phi_phi       <- array(0, dim = c(s, g, 1, 1))
  H_pi_pi         <- array(0, dim = c(s, g, 1, 1))

  fam_name <- family$family

  for (j in 1:s) {
    yj       <- Y[, j]
    beta0_j  <- beta0[j]
    lambda_j <- Lambda[j, ]
    phi_j    <- phi[j]
    a_phi    <- a_fun(fam_name, phi_j)

    for (k in 1:g) {

      beta_k <- B[k, ]

      ## η_{ijk}
      eta <- as.vector(beta0_j + X %*% beta_k + U %*% lambda_j)

      ## μ_{ijk}
      mu  <- family$linkinv(eta)

      ## μ'(η)
      mu_eta <- family$mu.eta(eta)

      ## V(μ)
      Vmu <- family$variance(mu)

      ## ∂²ℓ/∂η² = - V(μ) / (a(φ) * μ'(η)^2)
      w_eta <- - Vmu / (a_phi * mu_eta^2)

      ## β0β0
      H_beta0_beta0[j, k, 1, 1] <- sum(w_eta)

      ## BB
      H_B_B[j, k, , ] <- t(X * w_eta) %*% X

      ## ΛΛ
      H_Lambda_Lambda[j, k, , ] <- t(U * w_eta) %*% U

      ## β0B
      H_beta0_B[j, k, 1, ] <- colSums(X * w_eta)

      ## β0Λ
      H_beta0_Lambda[j, k, 1, ] <- colSums(U * w_eta)

      ## BΛ
      H_B_Lambda[j, k, , ] <- t(X * w_eta) %*% U

      ## φφ (zero for Poisson/Binomial)
      if (fam_name %in% c("poisson", "binomial")) {
        H_phi_phi[j, k, 1, 1] <- 0
      } else {
        stop("Hessian for φ_j not implemented for this family")
      }

      ## ππ
      H_pi_pi[j, k, 1, 1] <- -1 / (pi[k]^2)
    }
  }

  list(
    beta0_beta0   = H_beta0_beta0,
    B_B           = H_B_B,
    Lambda_Lambda = H_Lambda_Lambda,
    beta0_B       = H_beta0_B,
    beta0_Lambda  = H_beta0_Lambda,
    B_Lambda      = H_B_Lambda,
    phi_phi       = H_phi_phi,
    pi_pi         = H_pi_pi
  )
}

## helper: flatten (j,k) block scores into θ* vector
.flatten_score_jk <- function(j, k, score, s, g, p, d) {
  # θ* = (beta0[1:s], vec(B)[1:(g*p)], pi[1:g], vec(Lambda)[1:(s*d)], phi[1:s])
  P <- s + g * p + g + s * d + s
  v <- numeric(P)

  idx_beta0 <- seq_len(s)
  idx_B     <- (s + 1):(s + g * p)
  idx_pi    <- (s + g * p + 1):(s + g * p + g)
  idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
  idx_phi   <- (s + g * p + g + s * d + 1):(s + g * p + g + s * d + s)

  ## β0_j
  v[idx_beta0[j]] <- score$beta0[j, k, 1]

  ## β_k (length p), stored row-major in B
  B_offset <- (k - 1) * p
  v[idx_B[B_offset + seq_len(p)]] <- score$B[j, k, ]

  ## π_k
  v[idx_pi[k]] <- score$pi[j, k, 1]

  ## Λ_j (length d), stored row-major in Lambda
  Lam_offset <- (j - 1) * d
  v[idx_Lam[Lam_offset + seq_len(d)]] <- score$Lambda[j, k, ]

  ## φ_j
  v[idx_phi[j]] <- score$phi[j, k, 1]

  v
}

## helper: flatten (j,k) block Hessian into θ* × θ* matrix
.flatten_hessian_jk <- function(j, k, H, s, g, p, d) {
  P <- s + g * p + g + s * d + s
  M <- matrix(0, P, P)

  idx_beta0 <- seq_len(s)
  idx_B     <- (s + 1):(s + g * p)
  idx_pi    <- (s + g * p + 1):(s + g * p + g)
  idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
  idx_phi   <- (s + g * p + g + s * d + 1):(s + g * p + g + s * d + s)

  i_beta0 <- idx_beta0[j]
  i_pi    <- idx_pi[k]
  i_phi   <- idx_phi[j]

  i_B   <- idx_B[((k - 1) * p + 1):((k - 1) * p + p)]
  i_Lam <- idx_Lam[((j - 1) * d + 1):((j - 1) * d + d)]

  ## β0β0
  M[i_beta0, i_beta0] <- H$beta0_beta0[j, k, 1, 1]

  ## BB
  M[i_B, i_B] <- H$B_B[j, k, , ]

  ## ΛΛ
  M[i_Lam, i_Lam] <- H$Lambda_Lambda[j, k, , ]

  ## β0B and Bβ0
  M[i_beta0, i_B] <- H$beta0_B[j, k, 1, ]
  M[i_B, i_beta0] <- H$beta0_B[j, k, 1, ]

  ## β0Λ and Λβ0
  M[i_beta0, i_Lam] <- H$beta0_Lambda[j, k, 1, ]
  M[i_Lam, i_beta0] <- H$beta0_Lambda[j, k, 1, ]

  ## BΛ and ΛB
  M[i_B, i_Lam] <- H$B_Lambda[j, k, , ]
  M[i_Lam, i_B] <- t(H$B_Lambda[j, k, , ])

  ## φφ
  M[i_phi, i_phi] <- H$phi_phi[j, k, 1, 1]

  ## ππ
  M[i_pi, i_pi] <- H$pi_pi[j, k, 1, 1]

  M
}

#' Observed information matrix for a csam fit via Louis' identity
#'
#' θ* ordering: (beta0, vec(B), pi, vec(Lambda), phi)
#'
observed_information_csam <- function(object,
                                      score   = NULL,
                                      hessian = NULL) {

  if (!inherits(object, "csam"))
    stop("object must be a fitted csam model")

  Y      <- object$Y
  X      <- object$X
  U      <- object$U
  beta0  <- object$beta0
  B      <- object$B
  Lambda <- object$Lambda
  phi    <- object$phi
  pi     <- object$pi
  tau    <- object$tau

  s <- ncol(Y)
  g <- nrow(B)
  p <- ncol(X)
  d <- ncol(U)

  P <- s + g * p + g + s * d + s

  if (is.null(score))   score   <- score_all(object)
  if (is.null(hessian)) hessian <- hessian_all(object)

  I_obs <- matrix(0, P, P)

  ## Louis identity with H_{jk} = ∂²ℓ_{jk}/∂θ²:
  ##
  ## I_obs = -∑_{j,k} τ_{jk} H_{jk}
  ##         - ∑_j [ ∑_k τ_{jk} s_{jk}s_{jk}^T
  ##                 - (∑_k τ_{jk} s_{jk})(∑_k τ_{jk} s_{jk})^T ].

  for (j in 1:s) {

    Sj <- matrix(0, P, 1)  # ∑_k τ_{jk} s_{jk}
    Mj <- matrix(0, P, P)  # ∑_k τ_{jk} s_{jk}s_{jk}^T

    for (k in 1:g) {

      # scale tau (per-species average -> expected count)
      tau_jk <- tau[j, k] * n

      s_jk <- .flatten_score_jk(j, k, score, s, g, p, d)
      H_jk <- .flatten_hessian_jk(j, k, hessian, s, g, p, d)

      # -E_Z[∂²ℓ_c]
      I_obs <- I_obs - tau_jk * H_jk

      # accumulate for Var_Z term
      Sj <- Sj + tau_jk * s_jk
      Mj <- Mj + tau_jk * (s_jk %*% t(s_jk))
    }

    # subtract Var_Z(∂ℓ_c) contribution for species j
    I_obs <- I_obs - (Mj - Sj %*% t(Sj))
  }

  # symmetrise for numerical stability
  I_obs <- 0.5 * (I_obs + t(I_obs))

  I_obs
}

#' Covariance matrix for a csam fit via Louis' identity
#'
#' Returns vcov(θ*) where θ* = (beta0, vec(B), pi, vec(Lambda), phi).
#'
vcov_csam <- function(object,
                      score   = NULL,
                      hessian = NULL,
                      invert  = TRUE,
                      tol     = 1e-10,
                      marg_u = TRUE) {

  if (marg_u) {
    Iobs <- observed_information_csam_laplace_fixed(object, score, hessian)
  } else {
    Iobs <- observed_information_csam(object, score, hessian)
  }

  # symmetrise
  Iobs <- 0.5 * (Iobs + t(Iobs))

  if (!invert)
    return(Iobs)

  # safe inversion
  eig <- eigen(Iobs, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # regularise tiny eigenvalues
  vals[vals < tol] <- tol

  vc <- vecs %*% diag(1 / vals) %*% t(vecs)
  vc
}

#' Extract standard errors for θ* from a csam fit
#'
se_csam <- function(object, vcov_mat = NULL) {

  if (is.null(vcov_mat))
    vcov_mat <- vcov_csam(object)

  se <- sqrt(diag(vcov_mat))
  se
}

#' Reshape SE vector into csam-style parameter list (Option B)
#'
#' θ* ordering: (beta0, vec(B), pi, vec(Lambda), phi)
#' Returned list: (beta0, B, pi, U=NULL, Lambda, phi)
#'
se_list_csam <- function(object, se_vec = NULL) {

  if (is.null(se_vec))
    se_vec <- se_csam(object)

  Y      <- object$Y
  X      <- object$X
  U      <- object$U
  beta0  <- object$beta0
  B      <- object$B
  Lambda <- object$Lambda
  phi    <- object$phi
  pi     <- object$pi

  s <- ncol(Y)
  g <- nrow(B)
  p <- ncol(X)
  d <- ncol(U)

  # parameter indices
  idx_beta0 <- seq_len(s)
  idx_B     <- (s + 1):(s + g * p)
  idx_pi    <- (s + g * p + 1):(s + g * p + g)
  idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
  idx_phi   <- (s + g * p + g + s * d + 1):(s + g * p + g + s * d + s)

  # extract blocks
  se_beta0 <- se_vec[idx_beta0]

  se_B <- matrix(se_vec[idx_B], nrow = g, ncol = p, byrow = TRUE)

  se_pi <- se_vec[idx_pi]

  se_Lambda <- matrix(se_vec[idx_Lam], nrow = s, ncol = d, byrow = TRUE)

  se_phi <- se_vec[idx_phi]

  # Option B: include placeholder for U
  list(
    beta0  = se_beta0,
    B      = se_B,
    pi     = se_pi,
    U      = NULL,
    Lambda = se_Lambda,
    phi    = se_phi
  )
}

# Patched observed information with Laplace correction and diagnostics
observed_information_csam_laplace_fixed <- function(object,
                                                    score = NULL,
                                                    hessian = NULL,
                                                    tau = NULL,
                                                    include_trace = FALSE,
                                                    tol = 1e-10,
                                                    reg_eig = 1e-6,
                                                    verbose = TRUE) {
  if (!inherits(object, "csam")) stop("object must be a fitted csam model")

  # extract
  Y <- object$Y; X <- object$X; U <- object$U
  beta0 <- object$beta0; B <- object$B; Lambda <- object$Lambda
  phi <- object$phi; pi <- object$pi; fam <- object$family
  if (is.null(tau)) tau <- object$tau

  n <- nrow(Y); s <- ncol(Y); g <- nrow(B); p <- ncol(X); d <- ncol(U)
  P <- s + g * p + g + s * d + s

  if (is.null(score))   score   <- score_all(object)
  if (is.null(hessian)) hessian <- hessian_all(object)

  # ---- scale tau: per-species average -> expected counts per species-component
  tau_scaled <- tau * n

  # ---- helper flatteners (same as before) ----
  .flatten_score_jk_theta <- function(j, k, score, s, g, p, d) {
    P <- s + g * p + g + s * d + s
    v <- numeric(P)
    idx_beta0 <- seq_len(s)
    idx_B     <- (s + 1):(s + g * p)
    idx_pi    <- (s + g * p + 1):(s + g * p + g)
    idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
    idx_phi   <- (s + g * p + g + s * d + 1):(s + g * p + g + s * d + s)
    v[idx_beta0[j]] <- score$beta0[j, k, 1]
    B_offset <- (k - 1) * p
    v[idx_B[B_offset + seq_len(p)]] <- score$B[j, k, ]
    v[idx_pi[k]] <- score$pi[j, k, 1]
    Lam_offset <- (j - 1) * d
    v[idx_Lam[Lam_offset + seq_len(d)]] <- score$Lambda[j, k, ]
    v[idx_phi[j]] <- score$phi[j, k, 1]
    v
  }

  .flatten_hessian_jk_theta <- function(j, k, H, s, g, p, d) {
    P <- s + g * p + g + s * d + s
    M <- matrix(0, P, P)
    idx_beta0 <- seq_len(s)
    idx_B     <- (s + 1):(s + g * p)
    idx_pi    <- (s + g * p + 1):(s + g * p + g)
    idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
    idx_phi   <- (s + g * p + g + s * d + 1):(s + g * p + g + s * d + s)
    i_beta0 <- idx_beta0[j]
    i_B     <- idx_B[((k - 1) * p + 1):((k - 1) * p + p)]
    i_Lam   <- idx_Lam[((j - 1) * d + 1):((j - 1) * d + d)]
    i_phi   <- idx_phi[j]
    i_pi    <- idx_pi[k]
    M[i_beta0, i_beta0] <- H$beta0_beta0[j, k, 1, 1]
    M[i_B, i_B] <- H$B_B[j, k, , ]
    M[i_Lam, i_Lam] <- H$Lambda_Lambda[j, k, , ]
    M[i_beta0, i_B] <- H$beta0_B[j, k, 1, ]
    M[i_B, i_beta0] <- H$beta0_B[j, k, 1, ]
    M[i_beta0, i_Lam] <- H$beta0_Lambda[j, k, 1, ]
    M[i_Lam, i_beta0] <- H$beta0_Lambda[j, k, 1, ]
    M[i_B, i_Lam] <- H$B_Lambda[j, k, , ]
    M[i_Lam, i_B] <- t(H$B_Lambda[j, k, , ])
    M[i_phi, i_phi] <- H$phi_phi[j, k, 1, 1]
    M[i_pi, i_pi] <- H$pi_pi[j, k, 1, 1]
    M
  }

  # ---- compute scalars per (j,k) safely (uses your existing family analytic forms) ----
  .compute_scalars_per_jk <- function(object) {
    Y <- object$Y; X <- object$X; U <- object$U
    beta0 <- object$beta0; B <- object$B; Lambda <- object$Lambda
    phi <- object$phi; fam <- object$family
    n <- nrow(Y); s <- ncol(Y); g <- nrow(B)
    eta <- array(0, dim = c(n, s, g))
    mu  <- array(0, dim = c(n, s, g))
    mu_eta <- array(0, dim = c(n, s, g))
    delta <- array(0, dim = c(n, s, g))
    omega <- array(0, dim = c(n, s, g))
    kappa <- array(0, dim = c(n, s, g))
    fam_name <- fam$family
    for (k in 1:g) {
      beta_k <- B[k, ]
      for (j in 1:s) {
        eta_jk <- as.vector(beta0[j] + X %*% beta_k + U %*% Lambda[j, ])
        mu_jk  <- fam$linkinv(eta_jk)
        mu_eta_jk <- fam$mu.eta(eta_jk)
        eta[, j, k] <- eta_jk; mu[, j, k] <- mu_jk; mu_eta[, j, k] <- mu_eta_jk
        yj <- Y[, j]
        if (grepl("poisson", fam_name, ignore.case = TRUE)) {
          delta[, j, k] <- (yj - mu_jk) / mu_jk
          omega[, j, k] <- -1 / mu_jk
          kappa[, j, k] <- -(1 - 2 * mu_jk) / (mu_jk^3)
        } else if (grepl("binomial", fam_name, ignore.case = TRUE)) {
          delta[, j, k] <- (yj - mu_jk) / (mu_jk * (1 - mu_jk))
          omega[, j, k] <- -1 / (mu_jk * (1 - mu_jk))
          kappa[, j, k] <- - ( (1 - 2 * mu_jk) * (1 - 2 * mu_jk * (1 - mu_jk)) ) /
            (mu_jk^3 * (1 - mu_jk)^3)
        } else if (grepl("gaussian", fam_name, ignore.case = TRUE)) {
          a_phi <- phi[j]
          delta[, j, k] <- (yj - mu_jk) / a_phi
          omega[, j, k] <- -1 / a_phi
          kappa[, j, k] <- 0
        } else {
          eps <- 1e-6
          mu_eta2_jk <- (fam$mu.eta(eta_jk + eps) - fam$mu.eta(eta_jk - eps)) / (2 * eps)
          Vmu <- fam$variance(mu_jk)
          Vmu_p <- (fam$variance(mu_jk + eps) - fam$variance(mu_jk - eps)) / (2 * eps)
          a_phi <- ifelse(grepl("poisson|binomial", fam_name, ignore.case = TRUE), 1, phi[j])
          delta[, j, k] <- (yj - mu_jk) / (a_phi * mu_eta_jk)
          omega[, j, k] <- - Vmu / (a_phi * mu_eta_jk^2)
          kappa[, j, k] <- - (Vmu_p * mu_eta_jk - 2 * Vmu * mu_eta2_jk) / (a_phi * mu_eta_jk^4)
        }
      }
    }
    list(eta = eta, mu = mu, mu_eta = mu_eta,
         delta = delta, omega = omega, kappa = kappa)
  }

  # ---- assemble U-blocks and inverses using tau_scaled (expected counts) ----
  .assemble_U_blocks_tau_safe <- function(object, scalars, tau_scaled, reg_eig) {
    U <- object$U; Lambda <- object$Lambda
    n <- nrow(U); s <- ncol(object$Y); g <- nrow(object$B); d <- ncol(U)
    delta <- scalars$delta; omega <- scalars$omega
    grad_log_prior <- function(u_i) -u_i
    hess_log_prior  <- function(u_i) -diag(d)
    H_U_blocks <- vector("list", n); H_U_inv <- vector("list", n); grad_u <- vector("list", n)
    min_eig <- numeric(n)
    for (i in 1:n) {
      g_u <- grad_log_prior(U[i, ]); H_block <- hess_log_prior(U[i, ])
      for (j in 1:s) {
        lam_j <- Lambda[j, ]
        for (k in 1:g) {
          w <- tau_scaled[j, k]
          if (w == 0) next
          g_u <- g_u + w * delta[i, j, k] * lam_j
          H_block <- H_block + w * omega[i, j, k] * (lam_j %o% lam_j)
        }
      }
      H_block <- 0.5 * (H_block + t(H_block))
      eig <- eigen(H_block, symmetric = TRUE)
      vals <- eig$values; vecs <- eig$vectors
      min_eig[i] <- min(vals)
      vals[vals < reg_eig] <- reg_eig
      H_inv <- vecs %*% diag(1 / vals) %*% t(vecs)
      H_U_blocks[[i]] <- H_block; H_U_inv[[i]] <- H_inv; grad_u[[i]] <- g_u
    }
    list(H_U_blocks = H_U_blocks, H_U_inv = H_U_inv, grad_u = grad_u, min_eig = min_eig)
  }

  # ---- build cross matrices WITHOUT tau (w0 = 1) ----
  .build_cross_matrices_no_tau <- function(object, scalars) {
    Y <- object$Y; X <- object$X; U <- object$U; Lambda <- object$Lambda
    n <- nrow(Y); s <- ncol(Y); g <- nrow(object$B); d <- ncol(U); p <- ncol(X)
    P <- s + g * p + g + s * d + s
    delta <- scalars$delta; omega <- scalars$omega
    cross_site <- vector("list", n)
    for (i in 1:n) {
      cross_jk <- vector("list", s)
      for (j in 1:s) {
        cross_jk[[j]] <- vector("list", g)
        lam_j <- Lambda[j, ]
        for (k in 1:g) {
          C <- matrix(0, nrow = d, ncol = P)
          idx_beta0 <- seq_len(s)
          idx_B     <- (s + 1):(s + g * p)
          idx_Lam   <- (s + g * p + g + 1):(s + g * p + g + s * d)
          C[, idx_beta0[j]] <- C[, idx_beta0[j]] + omega[i, j, k] * lam_j
          B_offset <- (k - 1) * p
          C[, idx_B[B_offset + seq_len(p)]] <- C[, idx_B[B_offset + seq_len(p)]] +
            omega[i, j, k] * (lam_j %o% X[i, ])
          Lam_offset <- (j - 1) * d
          C[, idx_Lam[Lam_offset + seq_len(d)]] <- C[, idx_Lam[Lam_offset + seq_len(d)]] +
            ( delta[i, j, k] * diag(d) + omega[i, j, k] * (U[i, ] %o% lam_j) )
          cross_jk[[j]][[k]] <- C
        }
      }
      cross_site[[i]] <- cross_jk
    }
    cross_site
  }

  # ---- run computations ----
  scalars <- .compute_scalars_per_jk(object)
  Ublocks <- .assemble_U_blocks_tau_safe(object, scalars, tau_scaled, reg_eig)
  H_U_blocks <- Ublocks$H_U_blocks; H_U_inv <- Ublocks$H_U_inv; min_eig_HU <- Ublocks$min_eig
  cross_site <- .build_cross_matrices_no_tau(object, scalars)

  # precompute flattened complete-data scores and Hessians
  S_list <- vector("list", length = s * g); H_list <- vector("list", length = s * g)
  idx <- 1
  for (j in 1:s) for (k in 1:g) {
    S_list[[idx]] <- .flatten_score_jk_theta(j, k, score, s, g, p, d)
    H_list[[idx]] <- .flatten_hessian_jk_theta(j, k, hessian, s, g, p, d)
    idx <- idx + 1
  }

  # compute M_jk and diagnostics
  H_adj_list <- vector("list", length = s * g)
  M_norms <- matrix(0, nrow = s, ncol = g)
  H_norms <- matrix(0, nrow = s, ncol = g)
  idx <- 1
  for (j in 1:s) {
    for (k in 1:g) {
      Mjk <- matrix(0, P, P)
      for (i in 1:n) {
        Cik <- cross_site[[i]][[j]][[k]]
        if (!all(Cik == 0)) Mjk <- Mjk + t(Cik) %*% H_U_inv[[i]] %*% Cik
      }
      Hc <- H_list[[idx]]
      H_adj_list[[idx]] <- 0.5 * (Hc - Mjk + t(Hc - Mjk))
      M_norms[j, k] <- norm(Mjk, "F")
      H_norms[j, k] <- norm(Hc, "F")
      idx <- idx + 1
    }
  }

  # assemble observed information via Louis' identity (use tau_scaled here)
  I_obs <- matrix(0, P, P)
  idx <- 1
  for (j in 1:s) {
    Sj_sum <- matrix(0, P, 1); Mj <- matrix(0, P, P)
    for (k in 1:g) {
      tau_jk <- tau_scaled[j, k]
      s_jk <- S_list[[idx]]
      H_jk <- H_adj_list[[idx]]
      I_obs <- I_obs - tau_jk * H_jk
      Sj_sum <- Sj_sum + tau_jk * s_jk
      Mj <- Mj + tau_jk * (s_jk %*% t(s_jk))
      idx <- idx + 1
    }
    I_obs <- I_obs - (Mj - Sj_sum %*% t(Sj_sum))
  }

  I_obs <- 0.5 * (I_obs + t(I_obs))
  eigI <- eigen(I_obs, symmetric = TRUE)
  valsI <- eigI$values; vecsI <- eigI$vectors
  valsI[valsI < tol] <- tol
  I_obs_reg <- vecsI %*% diag(valsI) %*% t(vecsI)

  # diagnostics summary
  diag_list <- list(
    min_eig_HU = min_eig_HU,
    H_norms = H_norms,
    M_norms = M_norms,
    max_M_over_H = max(M_norms / (H_norms + 1e-12)),
    I_eigenvalues = valsI
  )

  if (verbose) {
    cat("Diagnostics (observed_information_csam_laplace_fixed):\n")
    cat("  min eigenvalues of H_U per site (first 10):", paste(round(min_eig_HU[1:min(10,length(min_eig_HU))], 3), collapse = ", "), "\n")
    cat("  max(M_jk / ||H_jk||_F) =", round(diag_list$max_M_over_H, 6), "\n")
    cat("  summary of I_obs eigenvalues: min =", round(min(valsI), 6), " max =", round(max(valsI), 6), "\n")
  }

  # list(I_obs = I_obs_reg, diagnostics = diag_list)
  return(I_obs_reg)
}
