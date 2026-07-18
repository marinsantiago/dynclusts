# ------------------------------------------------------------------------------
# Update tau2
# ------------------------------------------------------------------------------

update_tau2 <- function(gamma., K.check.inv, a.tau, b.tau) {
  quadform_gamma <- sum(gamma. * (K.check.inv %*% gamma.))
  a.post <- max(a.tau + length(gamma.) / 2, 1e-07)
  b.post <- max(b.tau + quadform_gamma / 2, 1e-07)
  max(1 / rgamma(1, shape = a.post, rate = b.post), 1e-10)
}
