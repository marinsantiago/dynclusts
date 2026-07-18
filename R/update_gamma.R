# ------------------------------------------------------------------------------
# Update gamma
# ------------------------------------------------------------------------------

update_gamma <- function(y, X.underbar,
                         theta, sigma2, beta., K.check.inv, tau2) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  # Linear predictor and residuals ---------------------------------------------
  linpred.const <- X.underbar %*% beta.
  dim(linpred.const) <- c(T.tot, n)
  y.check <- y - theta - linpred.const
  # Full conditional hyper-parameters ------------------------------------------
  K.inv <- (1 / tau2) * K.check.inv
  sigma2.inv <- 1 / sigma2
  sum.precisions <- diag(Rfast::colsums(sigma2.inv))
  post.precision <- K.inv + sum.precisions + (diag(n) * 1e-10)
  post.param.vector <- Rfast::colsums(sigma2.inv * y.check)
  # Sample from the full conditional -------------------------------------------
  gamma.out <- rand_mvnorm(post.precision, post.param.vector)
  gamma.out
}
