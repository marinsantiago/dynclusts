# ------------------------------------------------------------------------------
# Update psi
# ------------------------------------------------------------------------------

update_psi <- function(epsilon, alpha, lambda, psi,
                       sigma.psi_mh, iter, target.acc.rate = 0.44) {
  
  # Pre-compute constants ------------------------------------------------------
  dim.eps <- dim(epsilon)
  T.tot <- dim.eps[1]
  H1 <- dim.eps[2] # Implicitly, this is already H - 1
  T1 <- T.tot - 1
  H1T1 <- H1 * T1
  # Additional cache
  e.hat <- epsilon - (lambda / 2) * (1 - alpha)
  sum.e1T <- sum(Rfast::colsums((e.hat[c(1, T.tot),]^2) / lambda[c(1, T.tot),]))
  sum.e.mid <- Rfast::colsums(
    (e.hat[-c(1, T.tot),]^2) / lambda[-c(1, T.tot),]
  ) |> sum()
  sum.e.diff <- Rfast::colsums(
    (e.hat[-T.tot,] * e.hat[-1,]) / sqrt((lambda[-T.tot,] * lambda[-1,]))
  ) |> sum()
  # Evaluate log-full-conditional density and gradient at current psi ----------
  log_full_current <- logdens_psi(psi, H1T1, sum.e1T, sum.e.mid, sum.e.diff)
  lp.current <- log_full_current[1] # log-density
  gr.current <- log_full_current[2] # Gradient
  # Metropolis-Hastings algorithm with a Barker proposal -----------------------
  # 1. Sample z-noise
  z <- rnorm(1, 0, sd = sigma.psi_mh)
  # 2. Compute the probability of toggle
  p.toggle <- 1 / (1 + exp(-z * gr.current))
  # 3. Sample toggle
  toggle <- if (runif(1) < p.toggle) +1 else -1
  # 4. Sample the proposal 
  psi.star <- psi + toggle * z
  # 5. Compute M-H acceptance ratio
  if (abs(psi.star) < 1) { # Make sure the AR(1) process is stationary
    log_full_star <- logdens_psi(psi.star, H1T1, sum.e1T, sum.e.mid, sum.e.diff)
    lp.star <- log_full_star[1] # log-density
    gr.star <- log_full_star[2] # Gradient
    log.r <- lp.star - lp.current + as.list(softplus(
      c((psi - psi.star) * gr.current, (psi.star - psi) * gr.star)
    )) |> do.call("-", args = _)
  } else {
    log.r <- log(0) # Reject the proposal
  }
  # 6. Decision rule
  if (log(runif(1)) < log.r) psi <- psi.star
  # Robbins-Monro recursion ----------------------------------------------------
  rate <- min(exp(log.r), 1)
  sigma.psi <- exp(log(sigma.psi_mh) + iter^(-0.7) * (rate - target.acc.rate))
  list(psi = psi, sigma.psi_mh = sigma.psi)
}
