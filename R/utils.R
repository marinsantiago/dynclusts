# ------------------------------------------------------------------------------
# Utils and helpers
# ------------------------------------------------------------------------------

# Checks if an object inherits from "dynclusts".
is.dynclust <- \(object) inherits(object, "dynclusts")

# Checks if an object inherits from "dynclust.sim".
is.dynclust.sim <- \(object) inherits(object, "dynclust.sim")

# Softplus function to compute log(1 + exp(x))
softplus <- \(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Inverse logit
inv.logit <- \(x) 1 / (1 + exp(-x)) 

# Make a matrix  positive definite
make_posdef <- \(m) {
  d <- dim(m)[1]
  eigen_m <- eigen(m, symmetric = TRUE)
  evals_m <- eigen_m$values
  tol <- 2 * (d * max(abs(evals_m)) * .Machine$double.eps)
  delta = pmax(0, tol - evals_m)
  m + eigen_m$vectors %*% diag(delta, d) %*% t(eigen_m$vectors)
}

# Make a symmetric tri-diagonal matrix positive definite
make.posdef.tridiag <- \(diagonal, u.diagonal, buffer = 1e-6) {
  d <- diagonal
  e <- u.diagonal
  n <- length(d)
  d_new <- d
  for (i in seq_len(n)) {
    off_diag_sum <- 0
    if (i > 1) off_diag_sum <- off_diag_sum + abs(e[i - 1])
    if (i < n) off_diag_sum <- off_diag_sum + abs(e[i])
    if (d_new[i] <= off_diag_sum) d_new[i] <- off_diag_sum + buffer
  }
  d_new
}

# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------

# Check if the input is not a positive scalar
not.psc <- \(x) (!is.numeric(x)) || (length(x) != 1) || (x <= 0)

# Check if the input is not (1D) logical
not.logic <- \(x) (!is.logical(x)) || (length(x) != 1)

# Check if the input is not (1D) integer
not.int <- \(x) (x %% 1) != 0 || (x <= 0) || (length(x) != 1)

# Check if the input is not in the unit real line
not.unit <- \(x) (!is.numeric(x)) || (length(x) != 1) || (x <= 0) || (x >= 1)

# Check if the input is not a valid response matrix
not.y <- \(x) !is.numeric(x) || !is.matrix(x) || any(is.na(x))

# Check if the input is not a coordinate matrix: (Longitude, Latitude)
not.coord.mtrx <- \(x) {
  if (!is.matrix(x) || ncol(x) != 2L || !is.numeric(x) || any(is.na(x))) { 
    return(TRUE)
  }
  lon <- x[, 1]
  lat <- x[, 2]
  any(lat < -90 | lat > 90) || any(lon < -180 | lon > 180)
}

# Check if the input is not a valid list of matrices of covariates
not.X <- function(x, y) {
  if (!is.list(x)) return(TRUE)
  # All elements must be numeric matrices with no missing values
  if (!all(vapply(x, is.matrix, logical(1)))) return(TRUE)
  if (!all(vapply(x, is.numeric, logical(1)))) return(TRUE)
  if (any(vapply(x, anyNA, logical(1)))) return(TRUE)
  # Consistent dimensions
  ncols <- vapply(x, ncol, integer(1))
  nrows <- vapply(x, nrow, integer(1))
  if (length(unique(ncols)) != 1L || length(unique(nrows)) != 1L) return(TRUE)
  # Match rows with y
  nrows[1L] != nrow(y)
}

input.validation.dynclust <- \(y, coords, X, theta.0, sigma2.0, a.0, b.0, 
                               a.alpha, b.alpha, a.phi, b.phi, a.rho, b.rho, 
                               a.tau, b.tau, sigma.phi_mh, sigma.psi_mh, H, Sg,
                               post.pred, store.nu, store.epsilon, store.lambda, 
                               verbose, max.iters, burn.in, thin) {
  
  check <- !is.numeric(theta.0) || length(theta.0) != 1
  if (check) stop("theta.0 must be a real number")
  p.scalars <- c(
    "sigma2.0", "a.0", "b.0", "a.alpha", "b.alpha", "a.phi", "b.phi", "a.rho",
    "b.rho", "a.tau", "b.tau", "sigma.phi_mh", "sigma.psi_mh"
  )
  for (s in p.scalars) {
    if (not.psc(get(s))) stop(paste(s, " must be a positive scalar"))
  }
  logicals <- c(
    "Sg", "post.pred", "store.nu", "store.epsilon", "store.lambda", "verbose"
  )
  for (l in logicals) if (not.logic(get(l))) stop(paste(l, " must be logical"))
  integers <- c("H", "max.iters", "thin")
  for (i in integers) if (not.int(get(i))) stop(paste(i, " must be an integer"))
  check <- not.int(burn.in) || max.iters <= burn.in
  if (check) stop("burn.in must be a positive integer smaller than max.iters")
  rm(check, p.scalars, logicals, integers); gc()
  if (not.y(y)) stop("y must be a valid response matrix")
  if (not.X(X, y)) stop("X is not a valid list of covariates")
  if (not.coord.mtrx(coords)) stop("coords must be a valid coordinates matrix")
}

input.validation.dynclust.sim <- \(n.obs, time.points, n.clusts, n.covariates,
                                   country, tau2, rho2, varphi, sigma2, pi.stay, 
                                   balanced, spatial.effect, covariates.effect){
  
  for (s in c("tau2", "rho2", "varphi", "sigma2")) {
    if (not.psc(get(s))) stop(paste(s, " must be a positive scalar"))
  }
  integers <- c("n.obs", "time.points", "n.covariates", "n.clusts")
  for (i in integers) if (not.int(get(i))) stop(paste(i, " must be an integer"))
  logicals <- c("balanced", "spatial.effect", "covariates.effect")
  for (l in logicals) if (not.logic(get(l))) stop(paste(l, " must be logical"))
  if (not.unit(pi.stay)) stop("pi.stay must be in (0, 1)")
}
