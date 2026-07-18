# ------------------------------------------------------------------------------
# Update DP concentration parameter, alpha
# ------------------------------------------------------------------------------

# Stirling-gamma prior on alpha
update_alpha_sg <- function(K.tot, a.alpha, b.alpha, n) {
  # Full conditional hyper-parameters ------------------------------------------
  a.post <- a.alpha + K.tot
  b.post <- b.alpha + 1
  # Iterate over time points ---------------------------------------------------
  T.tot <- length(K.tot)
  alpha.out <- rep(NA, T.tot)
  for (t in seq_len(T.tot)) {
    alpha.out[t] <- ConjugateDP::rSg(1, a.post[t], b.post, n)
  }
  pmax(alpha.out, 1e-07)
}

# Gamma prior on alpha
update_alpha_g <- function(K.tot, a.alpha, b.alpha, n, alpha) {
  # Iterate over time points ---------------------------------------------------
  T.tot <- length(K.tot)
  alpha.out <- rep(NA, T.tot)
  for (t in seq_len(T.tot)) {
    # Auxiliary beta-distributed random variable
    eta <- rbeta(1, alpha[t] + 1, n)
    # Two-component mixture of gamma distributions
    pi.numerator <- a.alpha + K.tot[t] - 1
    pi.denominatpr <- n * (b.alpha - log(eta)) + pi.numerator
    pi_ <- pi.numerator / pi.denominatpr # Prob of the mixture components
    indicator <- rbinom(1, 1, pi_) # Indicator of mixture components
    alpha.out[t] <- if (indicator == 1) {
      # Sample from the first mixture component
      rgamma(1, shape = a.alpha + K.tot[t], rate = b.alpha - log(eta))
    } else {
      # Sample from the second mixture component
      rgamma(1, shape = a.alpha + K.tot[t] - 1, rate = b.alpha - log(eta))
    }
  }
  pmax(alpha.out, 1e-07)
}
