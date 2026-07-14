# ------------------------------------------------------------------------------
# Update rho2
# ------------------------------------------------------------------------------

update_rho2 <- function(beta., a.rho, b.rho) {
  a.post <- max(a.rho + length(beta.) / 2, 1e-07)
  b.post <- max(b.rho + sum(beta.^2) / 2, 1e-07)
  max(1 / rgamma(1, shape = a.post, rate = b.post), 1e-10)
}
