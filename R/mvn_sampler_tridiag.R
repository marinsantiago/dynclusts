# ------------------------------------------------------------------------------
# Multivariate normal sampler when precision matrix is tri-diagonal
# ------------------------------------------------------------------------------

# This function generates one random draw.
# param.vec: Vector such that Prec^{-1} %*% param.vec == mean.vector
# prec.diag: Main diagonal of the precision matrix as a vector
# prec.udiag: Upper diagonal of the precision matrix as a vector
rand_mvn_tridiag <- function(param.vec, prec.diag, prec.udiag) {
  
  # Pre-compute constants ------------------------------------------------------
  d <- length(param.vec)
  # Cholesky factorization of a tri-diagonal matrix ----------------------------
  cholesky.out <- mgcv::trichol(prec.diag, prec.udiag)
  L.diag <- cholesky.out$ld + 1e-08
  L.u.diag <- cholesky.out$sd
  # Get the mean vector via forward and backward solvers -----------------------
  intermediate.vec <- solve_lower_bidiagonal(L.diag, L.u.diag, param.vec)
  mean.vec <- solve_upper_bidiagonal(L.diag, L.u.diag, intermediate.vec)
  # Generate random draw from std normal ---------------------------------------
  z <- rnorm(d)
  # Solve system via backward solvers and return -------------------------------
  solve_upper_bidiagonal(L.diag, L.u.diag, z)
}
