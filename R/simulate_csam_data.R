simulate_csam_data <- function(n, s, g, d, p, true.pars = list(), covariate.matrix = NULL, family = poisson(), seed, plotting = TRUE) {

  set.seed(seed)

  # Default true parameter values
  beta0 <- if (!is.null(true.pars$beta0)) {true.pars$beta0} else {rnorm(s) / 4}
  B     <- if (!is.null(true.pars$B))     {true.pars$B}     else {B = matrix(rnorm(g * p), g, p); B[-1, ] = sweep(abs(B[-1, , drop = FALSE]), 2, -sign(B[1, ]), "*"); B}
  pi    <- if (!is.null(true.pars$pi))    {true.pars$pi}    else {rep(1/g, g)}
  U     <- if (!is.null(true.pars$U))     {true.pars$U}     else {matrix(rnorm(n * d), n, d)}
  Lambda<- if (!is.null(true.pars$Lambda)){true.pars$Lambda}else {Lambda = matrix(rnorm(s * d), s, d); Lambda[upper.tri(Lambda)] = 0; Lambda = Lambda %*% diag(sign(diag(Lambda)))}
  phi   <- if (!is.null(true.pars$phi))   {true.pars$phi}   else {rep(1, s)}

  # true indicators for species
  z0 = sample(1:g, s, replace = T, prob = pi)

  # env. covariate
  x <- if (!is.null(covariate.matrix)) {covariate.matrix} else {matrix(rnorm(n * p), n, p)}

  # arch. means
  eta_g = tcrossprod(x, B)

  # create a data frame for the sim
  dat = data.frame(sp = rep(1:s, each = n))
  dat$arch = z0[dat$sp] # arch. indicator
  dat$site = rep(1:n, s) # add in a site indicator
  dat$eta = NULL # species-specific mean
  dat$x = matrix(NA, n * s, p)
  dat$x = x[rep(seq_len(nrow(x)), s), ]
  for (j in 1:s) {
    dat$eta[dat$sp == j] = beta0[j] + eta_g[ , z0[j]] + (U %*% Lambda[j, ])
  }
  # obtain conditional mean based on family
  mu <- family$linkinv(dat$eta)

  # generate the response data based on the family
  if (family$family == "poisson") {

    dat$y <- rpois(nrow(dat), lambda = mu)

  } else if (family$family == "binomial") {

    dat$y <- rbinom(nrow(dat), size = 1, prob = mu)

  } else if (family$family == "gaussian") {

    dat$y <- rnorm(nrow(dat), mean = mu, sd = sqrt(dat$phi[dat$sp]))

  } else if (family$family == "Gamma") {

    shape <- 1 / dat$phi[dat$sp]
    scale <- mu / shape
    dat$y <- rgamma(nrow(dat), shape = shape, scale = scale)

  } else if (family$family == "inverse.gaussian") {

    if (!requireNamespace("statmod", quietly = TRUE)) {
      stop("statmod package required for inverse.gaussian simulation")
    }
    dat$y <- statmod::rinvgauss(nrow(dat), mean = mu,
                                dispersion = dat$phi[dat$sp])

  } else if (family$family == "nbinom2") {

    # NB2: Var = mu + phi * mu^2
    # MASS::rnegbin uses theta = 1/phi
    theta <- 1 / dat$phi[dat$sp]
    dat$y <- MASS::rnegbin(nrow(dat), mu = mu, theta = theta)

  } else {

    stop(sprintf("Simulation for family '%s' not implemented",
                 family$family))
  }

  # store the true parameters
  true.pars$beta0 = beta0
  true.pars$B = B
  true.pars$pi = pi
  true.pars$U = U
  true.pars$Lambda = Lambda
  true.pars$phi = phi

  attr(dat, "pars") = true.pars

  # Which looks something like this e.g.:
  if (plotting) {
    if (p > 1) {
      par(mfrow = c(1, p), mar = c(4.1, 4.1, 3.1, 0.5))
      for (k in 1:p) {
        with(dat, plot(x[,k], eta, col = arch, pch = 1, xlab = paste0("Env. Covariate ", k), ylab = expression(eta), main = ""))
        if (k == ceiling(p/2)) {
          legend(x = min(x[,k]) - (min(x[,k])*0.25), y = max(dat$eta) - (max(dat$eta)/8), legend = paste0("Archetype", 1:g), pch = 15, col = dat$arch, bty = "n", horiz = T, xpd = T)
        }
      }
      # legend("bottomleft", legend = c("Arch. 1", "Arch. 2"), pch = 1:2, bty = "n")
      # legend(x = -2.25, y = 6.5, legend = paste0("sp", 1:S), pch = 15, col = c25[1:S], bty = "n", horiz = T, xpd = T)
      par(mar = c(5.1, 4.1, 4.1, 2.1))
    } else {
      with(dat, plot(x, eta, col = arch, pch = 1, xlab = "Env. Covariate", ylab = expression(eta), main = ""))
      legend(x = -2.25, y = max(dat$eta), legend = paste0("Archetype", 1:g), pch = 15, col = 1:g, bty = "n", horiz = T, xpd = T)
    }
  }

  dat
}
