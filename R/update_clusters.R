# ------------------------------------------------------------------------------
# Update cluster indicators
# ------------------------------------------------------------------------------

update_clusters_disjoint <- function(y, X.underbar, H, omega, gamma.,
                                     beta., theta.candidates, sigma2) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  # Linear predictor  ----------------------------------------------------------
  linpred.const <- X.underbar %*% beta.
  dim(linpred.const) <- c(T.tot, n)
  linpred.const <- linpred.const + rep(gamma., each = T.tot) # Add theta later
  # Get the probabilities of each cluster --------------------------------------
  probs <- matrix(NA, nrow = n * T.tot, ncol = H)
  # Iterate over clusters
  for (k in seq_len(H)) {
    # Cluster-specific parameters
    current.mu <- linpred.const + theta.candidates[k] # Add theta
    probs[,k] <- rep(omega[,k], n) * as.vector(
      dnormal_disjoint(y, current.mu, sigma2)
    )
  }
  # Sample the cluster indicators s_it -----------------------------------------
  S.new <- int_sampling_rows(probs)
  matrix(S.new, nrow = T.tot, ncol = n, byrow = FALSE)
}

update_clusters_joint <- function(y, X.underbar, H, omega, gamma., beta.,
                                  theta.candidates, sigma2.candidates) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.y <- dim(y)
  T.tot <- dims.y[1]
  n <- dims.y[2]
  # Linear predictor  ----------------------------------------------------------
  linpred.const <- X.underbar %*% beta.
  dim(linpred.const) <- c(T.tot, n)
  linpred.const <- linpred.const + rep(gamma., each = T.tot) # Add theta later
  # Get the probabilities of each cluster --------------------------------------
  probs <- matrix(NA, nrow = n * T.tot, ncol = H)
  # Iterate over clusters
  for (k in seq_len(H)) {
    # Cluster-specific parameters
    current.mu <- linpred.const + theta.candidates[k] # Add theta
    probs[,k] <- rep(omega[,k], n) * as.vector(
      dnormal_joint(y, current.mu, sigma2.candidates[k])
    )
  }
  # Sample the cluster indicators s_it -----------------------------------------
  S.new <- int_sampling_rows(probs)
  matrix(S.new, nrow = T.tot, ncol = n, byrow = FALSE)
}
