# ------------------------------------------------------------------------------
# Update var-phi
# ------------------------------------------------------------------------------

update_varphi <- function(D, gamma., tau2, a.phi, b.phi, varphi, 
                          K.check, logdet.K.check, K.check.inv, 
                          sigma.phi_mh, iter, target.acc.rate = 0.44) {
  
  n <- length(gamma.)
  # Evaluate log-full-conditional density and gradient at current varphi -------
  log_full_current <- logdens_varphi(
    varphi, logdet.K.check, K.check.inv, K.check, D, gamma., tau2, a.phi, b.phi
  )
  lp.current <- log_full_current[1] # log-density
  gr.current <- log_full_current[2] # Gradient
  # Metropolis-Hastings algorithm with a Barker proposal -----------------------
  # 1. Sample z-noise
  z <- rnorm(1, 0, sd = sigma.phi_mh)
  # 2. Compute the probability of toggle
  p.toggle <- 1 / (1 + exp(-z * gr.current))
  # 3. Sample toggle
  toggle <- if (runif(1) < p.toggle) +1 else -1
  # 4. Sample the proposal 
  phi.star <- varphi + toggle * z
  # 5. Compute M-H acceptance ratio
  if (phi.star > 1e-07) { # Make sure the proposal is greater than zero
    # Update K.check matrix
    K.check.star <- exp(-D/(2 * phi.star^2)) + diag(n) * 1e-07
    K.check.star.chol <- tryCatch(
      { chol_float(K.check.star) }, error = \(e) NA
    )
    if (!any(is.na(K.check.star.chol))) { # Make sure the matrix is SPD
      logdet.K.check.star <- 2 * sum(log(diag(K.check.star.chol)))
      K.check.star.inv <- chol2inv(K.check.star.chol)
      log_full_star <- logdens_varphi(
        phi.star, logdet.K.check.star, K.check.star.inv, K.check.star,
        D, gamma., tau2, a.phi, b.phi
      )
      lp.star <- log_full_star[1] # log-density
      gr.star <- log_full_star[2] # Gradient
      log.r <- lp.star - lp.current + as.list(softplus(
        c((varphi - phi.star) * gr.current, (phi.star - varphi) * gr.star)
      )) |> do.call("-", args = _)
    } else {
      log.r <- log(0) # Reject the proposal
    }
  } else {
    log.r <- log(0) # Reject the proposal
  }
  # 6. Decision rule
  if (log(runif(1)) < log.r) {
    varphi <- phi.star # Update varphi
    K.check <- K.check.star # Update squared exponential kernel
    logdet.K.check <- logdet.K.check.star # Update log-determinant
    K.check.inv <- K.check.star.inv # Update inverse
  }
  # Robbins-Monro recursion ----------------------------------------------------
  rate <- min(exp(log.r), 1)
  sigma.phi <- exp(log(sigma.phi_mh) + iter^(-0.7) * (rate - target.acc.rate)) 
  list(
    varphi = varphi, K.check = K.check, logdet.K.check = logdet.K.check, 
    K.check.inv = K.check.inv, sigma.phi_mh = sigma.phi
  )
}
