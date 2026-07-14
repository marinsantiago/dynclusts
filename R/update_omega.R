# ------------------------------------------------------------------------------
# Update the stick-breaking weights, omega
# ------------------------------------------------------------------------------

update_omega <- function(S, alpha, lambda, psi, epsilon) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.S <- dim(S)
  T.tot <- dims.S[1]
  n <- dims.S[2]
  H1 <- ncol(epsilon) # H-1
  H <- H1 + 1
  psi2 <- psi^2 # Scalar
  Psi <- toeplitz(psi^(0:(T.tot - 1L)))
  # Iterate over clusters (up to H-1) ------------------------------------------
  epsilon.out <- lambda.out <- matrix(NA, nrow = T.tot, ncol = H1)
  for (k in seq_len(H1)) {
    # Get latent binomial variables --------------------------------------------
    I.idx <- S > k - 1
    Z.idx <- S == k
    m.idx <- Rfast::rowsums(I.idx)
    r.idx <- Rfast::rowsums(Z.idx)
    if (sum(m.idx) > 0) { # Make sure there is at least one time index
      # Step 1: Polya-gamma sampling -------------------------------------------
      xi <- rep(0, T.tot)
      # Two types of Polya-gamma samplers:
      # 1. Devroye-like (when m.k is a small integer)
      # 2. Saddle point approximation (otherwise)
      small.m <- (m.idx > 0 & m.idx <= 3); T.small <- sum(small.m)
      large.m <- (m.idx > 3L); T.large <- sum(large.m)
      if (T.small > 0) {
        xi[small.m] <- BayesLogit::rpg.devroye(
          T.small, m.idx[small.m], epsilon[small.m, k]
        )
      }
      if (T.large > 0) {
        xi[large.m] <- BayesLogit::rpg.sp(
          T.large, m.idx[large.m], epsilon[large.m, k]
        )
      }
      xi <- pmax(xi, 1e-08)
      # Step 2: Metropolis-Hastings step  --------------------------------------
      idx.omit <- m.idx == 0
      idx.keep <- !idx.omit
      T.k <- sum(idx.keep)
      xi.k <- xi[idx.keep]
      m.k <- m.idx[idx.keep]
      r.k <- r.idx[idx.keep]
      Psi.k <- Psi[idx.keep, idx.keep] |> as.matrix()
      alpha.k <- alpha[idx.keep]
      lambda.k <- lambda[idx.keep, k]
      # Generate the proposals
      proposals <- rand_polya_proposal(rep(1, T.k), alpha.k, lambda.k)
      lambda.star.k <- pmax(proposals[1,], 1e-08)
      # Evaluate "L" densities
      adj <- diag(T.k) * 1e-07
      xi.inv <- 1 / xi.k
      eval.point <- xi.inv * (r.k - m.k / 2)
      center.current <- eval.point - (lambda.k * (1 - alpha.k) / 2)
      center.star <- eval.point - (lambda.star.k * (1 - alpha.k) / 2)
      xi.inv.mat <- if (T.k > 1) diag(xi.inv) else xi.inv
      cov.current <- Psi.k * tcrossprod(sqrt(lambda.k)) + xi.inv.mat + adj
      cov.star <- Psi.k * tcrossprod(sqrt(lambda.star.k)) + xi.inv.mat + adj
      log.L.current <- lambda_ratio_log_l(center.current, cov.current)
      log.L.star <- lambda_ratio_log_l(center.star, cov.star)
      # M-H acceptance ratio
      a.prime <- proposals[2,]
      b.prime <- proposals[3,]
      log.L.diff <- log.L.star - log.L.current
      oo. <- sum((alpha.k - a.prime * b.prime) * (lambda.k - lambda.star.k)) / 2
      log.r <- log.L.diff + oo.
      # Decision rule
      lambda.out[idx.keep, k] <- if (log(runif(1)) < log.r) {
        lambda.star.k
      } else { lambda.k }
      # Sample the remaining lambdas, if any, straight from the prior
      T.omit <- sum(idx.omit)
      if (T.omit > 0) {
        lambda.out[idx.omit, k] <- vapply(
          seq_len(T.omit), \(t) rpolya(1, alpha[t]), 1
        )
      }
      # Step 3: MVN sampler ----------------------------------------------------
      # Note that Psi is an AR(1) matrix, so its inverse is tri-diagonal
      # Main diagonal:
      Psi.inv.diag <- rep(NA, T.tot) 
      Psi.inv.diag[c(1, T.tot)] <- 1 / (1 - psi2)
      Psi.inv.diag[-c(1, T.tot)] <- (1 + psi2) / (1 - psi2)
      # Super diagonal:
      Psi.inv.u.diag <- rep(-psi, T.tot - 1) / (1 - psi2)
      # Compute Lambda^{-1/2} %*% Psi^{-1} %*% Lambda^{-1/2}, which is tri-diag
      lambda.sqrt.inv <- 1 / sqrt(pmax(lambda.out[,k], 1e-07))
      m.diag <- lambda.sqrt.inv^2 * Psi.inv.diag 
      m.u.diag <- lambda.sqrt.inv[-T.tot] * Psi.inv.u.diag * lambda.sqrt.inv[-1]
      # Since Xi is diagonal, the precision Xi + M is also tri-diagonal
      # For the diagonal entries in which m.k != 0, add Xi (from step 1)
      xi.diag <- rep(0, T.tot)
      xi.diag[idx.keep] <- xi.k
      post.prec.diag <- m.diag + xi.diag
      post.prec.u.diag <- m.u.diag
      # Compute the parameter vector
      alpha.diff <- (1 - alpha)
      P.inv.alph <- Psi.inv.diag * alpha.diff # Psi^{-1} %*% (1 - alpha)
      P.inv.alph[-1] <- P.inv.alph[-1] + Psi.inv.u.diag * alpha.diff[-T.tot]
      P.inv.alph[-T.tot] <- P.inv.alph[-T.tot] + Psi.inv.u.diag * alpha.diff[-1]
      binomial.counts <- rep(0, T.tot) # For the entries in which m.k != 0
      binomial.counts[idx.keep] <- r.k - m.k / 2
      param.vec <- binomial.counts + P.inv.alph / 2
      # Sample from the MVN distribution exploiting the tri-diagonal precision
      eps.new <- tryCatch(
        { rand_mvn_tridiag(param.vec, post.prec.diag, post.prec.u.diag) }, 
        error = \(e) {
          # If needed, enforce SPD precision matrix
          rand_mvn_tridiag(
            param.vec, 
            make.posdef.tridiag(post.prec.diag, post.prec.u.diag), 
            post.prec.u.diag
          )
        }
      )
      epsilon.out[,k] <- eps.new
    } else {
      # If all the indices are empty, sample epsilon_k straight from the prior
      lambda.m.zero <- vapply(T.tot, \(t) rpolya(1, alpha[t]), 1) # prior draw 
      lambda.out[, k] <- lambda.m.zero
      # Note that Psi is an AR(1) matrix, so its inverse is tri-diagonal
      # Main diagonal:
      Psi.inv.diag <- rep(NA, T.tot) 
      Psi.inv.diag[c(1, T.tot)] <- 1 / (1 - psi2)
      Psi.inv.diag[-c(1, T.tot)] <- (1 + psi2) / (1 - psi2)
      # Super diagonal:
      Psi.inv.u.diag <- rep(-psi, T.tot - 1) / (1 - psi2)
      # Again, note that the prior precision is tri-diagonal!
      la.sqrt.inv <- 1 / sqrt(pmax(lambda.out[,k], 1e-07))
      prec.diag <- la.sqrt.inv^2 * Psi.inv.diag + 1e-06
      prec.u.diag <- la.sqrt.inv[-T.tot] * Psi.inv.u.diag * la.sqrt.inv[-1]
      # Compute the parameter vector
      alpha.diff <- (1 - alpha)
      P.inv.alph <- Psi.inv.diag * alpha.diff # Psi^{-1} %*% (1 - alpha)
      P.inv.alph[-1] <- P.inv.alph[-1] + Psi.inv.u.diag * alpha.diff[-T.tot]
      P.inv.alph[-T.tot] <- P.inv.alph[-T.tot] + Psi.inv.u.diag * alpha.diff[-1]
      param.vec <- P.inv.alph / 2
      # Sample from the MVN distribution exploiting the tri-diagonal precision
      epsilon.out[,k] <- rand_mvn_tridiag(param.vec, prec.diag, prec.u.diag)
    }
  }
  # Compute nu -----------------------------------------------------------------
  nu.out <- matrix(NA, nrow = T.tot, ncol = H)
  nu.out[,seq_len(H1)] <- inv.logit(epsilon.out)
  nu.out[,H] <- 1 # nu for the H-th cluster is always one
  # Compute omega --------------------------------------------------------------
  omega.out <- matrix(NA, nrow = T.tot, ncol = H)
  omega.out[,1] <- nu.out[,1]
  for (t in seq_len(T.tot)) {
    for (k in 2:H) {
      omega.out[t,k] <- nu.out[t,k] * prod(1 - nu.out[t,1:(k-1)])
    }
  }
  # Prepare the returns --------------------------------------------------------
  list(
    omega.out = omega.out, nu.out = nu.out, 
    epsilon.out = epsilon.out, lamda.out = lambda.out
  )
}
