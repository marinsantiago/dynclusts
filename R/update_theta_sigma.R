# ------------------------------------------------------------------------------
# Update theta and sigma2
# ------------------------------------------------------------------------------

update_theta_sigma2 <- function(y, X.underbar, S, H, beta., gamma., 
                                theta, sigma2, theta.0, sigma2.0, a.0, b.0) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  sigma.0 <- sqrt(sigma2.0)
  i.sigma2.0 <- 1 / sigma2.0
  theta.0.adj <- theta.0 * i.sigma2.0
  adj <- 1e-25
  # Linear predictor and residuals ---------------------------------------------
  # linpred.underbar <- X.underbar %*% beta. + rep(gamma., each = T.tot)
  # y.tilde <- y - matrix(linpred.underbar, ncol = n, byrow = FALSE)
  linpred.underbar <- drop(X.underbar %*% beta.)
  dim(linpred.underbar) <- c(T.tot, n)
  y.tilde <- y - linpred.underbar - rep(gamma., each = T.tot) 
  # Sample from the full conditionals of theta and sigma2 ----------------------
  theta.out <- sigma2.out <- rep(NA, H)
  # Iterate over clusters
  for (k in seq_len(H)) {
    y.clust <- y.tilde[S == k]
    n.clust <- length(y.clust)
    if (n.clust == 0) {
      # Sample from the prior centering measure
      theta.out[k] <- rnorm(1, mean = theta.0, sd = sigma.0)
      sigma2.out[k] <- max(1 / rgamma(1, shape = a.0, rate = b.0), adj)
    } else {
      # Update theta conditional on sigma2
      current.sigma2 <- sigma2[S == k][1]
      i.current.sigma2 <- 1 / current.sigma2
      sigma2.hat <- max(1 / (i.sigma2.0 + n.clust * i.current.sigma2), adj)
      theta.hat <- sigma2.hat * (theta.0.adj + i.current.sigma2 * sum(y.clust))
      curr.theta <- rnorm(1, mean = theta.hat, sd = sqrt(sigma2.hat))
      theta.out[k] <- curr.theta
      # Update sigma2 conditional on theta
      a.hat <- max(a.0 + n.clust / 2, adj)
      b.hat <- max(b.0 + sum((y.clust - curr.theta)^2) / 2, adj)
      sigma2.out[k] <- max(1 / rgamma(1, shape = a.hat, rate = b.hat), adj)
    }
  }
  list(theta.out = theta.out, sigma2.out = sigma2.out)
}
