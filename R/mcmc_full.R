# ------------------------------------------------------------------------------
# MCMC algorithm
# ------------------------------------------------------------------------------

mcmc_sampler <- function(y, coords, X, theta.0, sigma2.0, a.0, b.0, 
                         a.alpha, b.alpha, a.phi, b.phi, a.rho, b.rho, 
                         a.tau, b.tau, sigma.phi_mh, sigma.psi_mh, H, Sg, 
                         post.pred, store.nu, store.epsilon, store.lambda, 
                         verbose, max.iters, burn.in, thin) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  adj <- 1e-06
  # Covariates
  X.underbar <- do.call(rbind, X)
  p <- ncol(X.underbar)
  # Squared distances
  D <- if (!is.null(coords)) (haver_dist(coords) / 1000)^2 else matrix(0, n, n)
  # Initialize data structures to store posterior draws ------------------------
  # Number of posterior draws after burn-in and thinning
  dd <- floor((max.iters - burn.in) / thin)
  psi.out <- rho2.out <- varphi.out <- tau2.out <- rep(NA, dd)
  alpha.out <- K.tot.out <- matrix(nrow = dd, ncol = T.tot)
  beta.out <- matrix(nrow = dd, ncol = p)
  gamma.out <- matrix(nrow = dd, ncol = n)
  omega.out <- array(dim = c(T.tot, H, dd))
  theta.out <- sigma2.out <- S.out <- loglik.out <- array(dim = c(T.tot, n, dd))
  if (store.nu) nu.out <- array(dim = c(T.tot, H, dd))
  if (store.epsilon) epsilon.out <- array(dim = c(T.tot, H - 1, dd))
  if (store.lambda) lambda.out <- array(dim = c(T.tot, H - 1, dd))
  if (post.pred) pred.out <- array(dim = c(T.tot, n, dd))
  if (post.pred) X.time <- aperm(simplify2array(X), c(3, 2, 1)) # n x p x T.tot
  current.save <- 0
  # Initialize model parameters ------------------------------------------------
  psi <- runif(1)
  rho2 <- max(1 / rgamma(1, shape = a.rho, rate = b.rho), adj)
  beta. <- rnorm(p)
  gamma. <- rnorm(n)
  tau2 <- max(1 / rgamma(1, shape = a.tau, rate = b.tau), adj)
  varphi <- max(rgamma(1, shape = a.phi, rate = b.phi), adj)
  alpha <- pmax(ConjugateDP::rSg(T.tot, a.alpha, b.alpha, n), adj)
  # Squared exponential kernel
  K.check <- make_posdef(exp(-D / (2 * varphi^2)) + diag(n) * adj)
  K.check.chol <- chol(K.check)
  logdet.K.check <- 2 * sum(log(diag(K.check.chol)))
  K.check.inv <- chol2inv(K.check.chol)
  # Autoregressive stick-breaking weights
  epsilon <- lambda <- matrix(nrow = T.tot, ncol = H - 1)
  for (k in seq_len(H - 1)) {
    epsilon[,k] <- arima.sim(n = T.tot, model = list(ar = psi))
    lambda[,k] <- vapply(alpha, \(a) rpolya(1, a), 1)
  }
  nu <- matrix(nrow = T.tot, ncol = H)
  nu[,seq_len(H - 1)] <- inv.logit(epsilon)
  nu[,H] <- 1 # nu for the H-th cluster is always one
  omega <- matrix(nrow = T.tot, ncol = H)
  omega[,1] <- nu[,1]
  for (t in seq_len(T.tot)) {
    for (k in 2:H) {
      omega[t, k] <- nu[t, k] * prod(1 - nu[t, 1:(k-1)])
    }
  }
  # Cluster allocations
  S <- matrix(nrow = T.tot, ncol = n)
  for (t in seq_len(T.tot)) {
    S[t,] <- sample.int(H, n, replace = T, prob = omega[t,])
  }
  K.tot <- row_unique_counts(S)
  # Generate theta and sigma2
  theta.candidates <- seq(min(y) + adj, max(y) - adj, length.out = H)
  sigma2.candidates <- pmax(1 / rgamma(H, shape = a.0, rate = b.0), adj)
  theta <- sigma2 <- matrix(nrow = T.tot, ncol = n)
  for (i in seq_len(n)) {
    theta[,i] <- theta.candidates[S[,i]]
    sigma2[,i] <- sigma2.candidates[S[,i]]
  }
  # MCMC algorithm -------------------------------------------------------------
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = " dynclusts MCMC algorithm: [:bar] :percent in :elapsed", 
      total = max.iters, clear = FALSE
    )
  }
  start <- Sys.time() # Start wall-clock time
  for (iter in seq_len(max.iters)) {
    if (verbose) pb$tick()
    # Update cluster indicators
    S <- update_clusters(
      y = y, X.underbar = X.underbar, H = H, omega = omega, 
      gamma. = gamma., beta. = beta., theta.candidates, sigma2.candidates
    )
    K.tot <- row_unique_counts(S)
    # Update theta and sigma2
    theta.sigma2.candidates <- update_theta_sigma2(
      y = y, X.underbar = X.underbar, S = S, H = H, beta. = beta., 
      gamma. = gamma., theta = theta, sigma2 = sigma2, theta.0 = theta.0,
      sigma2.0 = sigma2.0, a.0 = a.0, b.0 = b.0
    )
    theta.candidates <- theta.sigma2.candidates$theta.out
    sigma2.candidates <- theta.sigma2.candidates$sigma2.out
    for (i in seq_len(n)) {
      theta[,i] <- theta.candidates[S[,i]]
      sigma2[,i] <- sigma2.candidates[S[,i]]
    }
    # Update beta
    beta. <- update_beta(
      y = y, X.underbar = X.underbar, theta = theta, 
      sigma2 = sigma2, gamma. = gamma., rho2 = rho2
    )
    # Update gamma
    gamma. <- update_gamma(
      y = y, X.underbar = X.underbar, theta = theta, sigma2 = sigma2,
      beta. = beta., K.check.inv = K.check.inv, tau2 = tau2
    )
    # Update stick-breaking weights
    stick.weights <- update_omega(
      S = S, alpha = alpha, lambda = lambda, psi = psi, epsilon = epsilon
    )
    nu <- stick.weights$nu.out
    omega <- stick.weights$omega.out
    lambda <- stick.weights$lambda.out
    epsilon <- stick.weights$epsilon.out
    # Update rho2
    rho2 <- update_rho2(beta. = beta., a.rho = a.rho, b.rho = b.rho)
    # Update tau2
    tau2 <- update_tau2(
      gamma. = gamma., K.check.inv = K.check.inv, a.tau = a.tau, b.tau = b.tau
    )
    # Update var-phi
    varphi.iter <- update_varphi(
      D = D, gamma. = gamma., tau2 = tau2, a.phi = a.phi, b.phi = b.phi,
      varphi = varphi, K.check = K.check, logdet.K.check = logdet.K.check,
      K.check.inv = K.check.inv, sigma.phi_mh = sigma.phi_mh, iter = iter
    ) 
    varphi <- varphi.iter$varphi
    K.check <- varphi.iter$K.check
    K.check.inv <- varphi.iter$K.check.inv
    logdet.K.check <- varphi.iter$logdet.K.check
    sigma.phi_mh <- varphi.iter$sigma.phi_mh
    # Update alpha
    alpha <- if (Sg) {
      # Stirling-gamma prior
      update_alpha_sg(
        K.tot = K.tot, a.alpha = a.alpha, b.alpha = b.alpha, n = n
      )
    } else {
      # Gamma prior
      update_alpha_g(
        K.tot = K.tot, a.alpha = a.alpha, 
        b.alpha = b.alpha, n = n, alpha = alpha
      )
    }
    # Update psi
    psi.iter <- update_psi(
      epsilon = epsilon, alpha = alpha, lambda = lambda, 
      psi = psi, sigma.psi_mh = sigma.psi_mh, iter = iter
    )
    psi <- psi.iter$psi
    sigma.psi_mh <- psi.iter$sigma.psi_mh
    # Store posterior draws after burn-in and thinning
    if (iter > burn.in) {
      if (iter %% thin == 0) {
        current.save <- current.save + 1L
        # Current log-likelihood and posterior predictive draws
        current.llike <- matrix(nrow = T.tot, ncol = n)
        if (post.pred) current.pred <- matrix(nrow = T.tot, ncol = n)
        for (t in seq_len(T.tot)) {
          mu <- theta[t,] + gamma. + X.time[,,t] %*% beta.
          sigma <- sqrt(sigma2[t,]) 
          current.llike[t,] <- dnorm(y[t,], mean = mu, sd = sigma, log = TRUE)
          if (post.pred) current.pred[t,] <- rnorm(n, mean = mu, sd = sigma)
        }
        loglik.out[,,current.save] <- current.llike
        if (post.pred) pred.out[,,current.save] <- current.pred
        if (store.nu) nu.out[,,current.save] <- nu
        if (store.epsilon) epsilon.out[,,current.save] <- epsilon
        if (store.lambda) lambda.out[,,current.save] <- lambda
        sigma2.out[,,current.save] <- sigma2
        omega.out[,,current.save] <- omega
        theta.out[,,current.save] <- theta
        gamma.out[current.save, ] <- gamma.
        K.tot.out[current.save, ] <- K.tot
        alpha.out[current.save, ] <- alpha
        beta.out[current.save, ] <- beta.
        varphi.out[current.save] <- varphi
        rho2.out[current.save] <- rho2
        tau2.out[current.save] <- tau2
        psi.out[current.save] <- psi
        S.out[,,current.save] <- S
      }
    }
    if (iter %% 5000 == 0) gc() # Collect cache
  }
  end <- Sys.time()
  elapsed <- difftime(end, start, units = "mins")
  # Returns --------------------------------------------------------------------
  mcmc.out <- list(
    theta.post = theta.out, sigma2.post = sigma2.out, omega.post = omega.out,
    psi.post = psi.out, beta.post = beta.out, rho2.post = rho2.out, 
    K.post = K.tot.out, tau2.post = tau2.out, alpha.post = alpha.out,
    gamma.post = gamma.out, varphi.post = varphi.out, S.post = S.out,
    loglik = loglik.out, elapsed = elapsed
  )
  if (store.nu) mcmc.out[["nu.post"]] <- nu.out
  if (store.epsilon) mcmc.out[["epsilon.post"]] <- epsilon.out
  if (store.lambda) mcmc.out[["lambda.post"]] <- lambda.out
  if (post.pred) mcmc.out[["post.pred"]] <- pred.out
  gc()
  mcmc.out
}
