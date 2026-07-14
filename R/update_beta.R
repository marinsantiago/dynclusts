# ------------------------------------------------------------------------------
# Update beta
# ------------------------------------------------------------------------------

update_beta <- function(y, X.underbar, theta, sigma2, gamma., rho2) {
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  p <- ncol(X.underbar)
  adj <- 1e-10
  # Linear predictor and residuals
  # linpred <- theta + matrix(gamma., nrow = T.tot, ncol = n, byrow = TRUE)
  # y.underbar <- c(y - linpred)
  y.underbar <- as.vector(y - theta) - rep(gamma., each = T.tot)
  # Full conditional hyper-parameters ------------------------------------------
  Sigma.underbar.inv <- 1 / as.vector(sigma2) # Just store the diagonal elements
  Xt_Sigma.inv <- Rfast::transpose(X.underbar * Sigma.underbar.inv)
  post.precision <- diag(p) * (1 + adj) / rho2 + Xt_Sigma.inv %*% X.underbar
  post.param.vector <- Xt_Sigma.inv %*% y.underbar
  # Sample from the full conditional distribution ------------------------------
  beta.out <- rand_mvnorm(post.precision, post.param.vector)
  beta.out
}
