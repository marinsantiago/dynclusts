# ------------------------------------------------------------------------------
# Update the stick-breaking weights, i.e., omega
# ------------------------------------------------------------------------------

update_omega <- function(S, alpha, lambda, psi, epsilon) {
  
  # Pre-compute constants ------------------------------------------------------
  dims.S <- dim(S)
  T.tot <- dims.S[1]
  n <- dims.S[2]
  H1 <- ncol(epsilon) # H-1
  H <- H1 + 1
  psi2 <- psi^2 # scalar
  Psi <- toeplitz(psi^(0:(T.tot - 1L)))
  # Iterate over clusters (up to H-1) ------------------------------------------
  epsilon.out <- lambda.out <- matrix(NA, nrow = T.tot, ncol = H1)
  for (k in seq_len(H1)) {
    # Get latent binomial variables --------------------------------------------
    I.idx <- S > k - 1
    Z.idx <- S == k
    m_k <- Rfast::rowsums(I.idx)
    r_k <- Rfast::rowsums(Z.idx)
    # Step 1: Polya-gamma sampling ---------------------------------------------
    xi <- rep(0, T.tot) # Sets to zero entries for which m_k == 0
    # Two types of Polya-gamma samplers:
    # (a). Devroye-like (when m.k is a small integer)
    # (b). Saddle point approximation (otherwise)
    small.m <- (m_k > 0 & m_k <= 3); T.small <- sum(small.m)
    large.m <- (m_k > 3L); T.large <- sum(large.m)
    if (T.small > 0) {
      xi[small.m] <- pmax(
        BayesLogit::rpg.devroye(T.small, m_k[small.m], epsilon[small.m, k]), 
        1e-20
      )
    }
    if (T.large > 0) {
      xi[large.m] <- pmax(
        BayesLogit::rpg.sp(T.large, m_k[large.m], epsilon[large.m, k]), 1e-20
      )
    }
    # Step 2: Metropolis-Hastings step  ----------------------------------------
    idx.omit <- m_k == 0; T.omit <- sum(idx.omit)
    idx.keep <- !idx.omit; T.k <- sum(idx.keep)
    if (T.omit > 0) {
      # Sample the lambdas for which m_k == 0 straight from the prior
      alpha.omit <- alpha[idx.omit]
      lambda.out[idx.omit, k] <- vapply(alpha.omit, \(a) rpolya(1, a), 1)
    }
    if (T.k > 0) {
      # M-H update for lambdas
      xi.keep <- xi[idx.keep]
      m.keep <- m_k[idx.keep]
      r.keep <- r_k[idx.keep]
      alpha.keep <- alpha[idx.keep]
      lambda.keep <- lambda[idx.keep, k]
      Psi.keep <- as.matrix(Psi[idx.keep, idx.keep])
      # Generate the proposals
      proposals <- rand_polya_proposal(rep(1, T.k), alpha.keep, lambda.keep)
      lambda.star.k <- pmax(proposals[1,], 1e-10)
      a.prime <- proposals[2,]
      b.prime <- proposals[3,]
      # Evaluate "L" functions
      adj <- diag(T.k) * 1e-10
      xi.inv <- 1 / xi.keep
      eval.point <- xi.inv * (r.keep - m.keep / 2)
      center.current <- eval.point - (lambda.keep * (1 - alpha.keep) / 2)
      center.star <- eval.point - (lambda.star.k * (1 - alpha.keep) / 2)
      xi.inv.mat <- if (T.k > 1) diag(xi.inv) else xi.inv
      cov.current <- Psi.keep * tcrossprod(sqrt(lambda.keep)) + xi.inv.mat + adj
      cov.star <- Psi.keep * tcrossprod(sqrt(lambda.star.k)) + xi.inv.mat + adj
      L.current <- lambda_ratio_log_l(center.current, cov.current)
      L.star <- lambda_ratio_log_l(center.star, cov.star)
      # M-H acceptance ratio
      L.diff <- L.star - L.current
      temp.. <- sum(
        (alpha.keep - a.prime * b.prime) * (lambda.keep - lambda.star.k) / 2
      )
      log.r <- L.diff + temp..
      # Decision rule
      lambda.MH <- if (log(runif(1)) < log.r) lambda.star.k else lambda.keep
      lambda.out[idx.keep, k] <- lambda.MH
    }
    # Step 3: MVN sampler ------------------------------------------------------
    # Note that Psi is a toeplitz matrix, so its inverse is tri-diagonal
    # Main diagonal:
    Psi.inv.diag <- rep(NA, T.tot)
    Psi.inv.diag[c(1, T.tot)] <- 1 / (1 - psi2)
    Psi.inv.diag[-c(1, T.tot)] <- (1 + psi2) / (1 - psi2)
    # Super diagonal:
    Psi.inv.u.diag <- rep(-psi, T.tot - 1) / (1 - psi2)
    # Compute Lambda^{-1/2} %*% Psi^{-1} %*% Lambda^{-1/2}, which is tri-diag
    lambda.sqrt.inv <- 1 / sqrt(pmax(lambda.out[,k], 1e-07))
    mm.diag <- lambda.sqrt.inv^2 * Psi.inv.diag 
    mm.u.diag <- lambda.sqrt.inv[-1] * Psi.inv.u.diag * lambda.sqrt.inv[-T.tot]
    # Since Xi is diagonal, the precision Xi + M is also tri-diagonal
    post.prec.diag <- mm.diag + xi
    post.prec.u.diag <- mm.u.diag
    # Compute the parameter vector
    bin.counts <- r_k - m_k / 2
    l.alpha <- lambda.out[,k] * (1 - alpha) / 2 # diag(lambda) %*% (1 - alpha)/2
    # Compute (Lambda^{-1/2} %*% Psi^{-1} %*% Lambda^{-1/2}) %*% l.alpha
    l.alpha.adj <- l.alpha * mm.diag
    l.alpha.adj[-T.tot] <- l.alpha.adj[-T.tot] + mm.u.diag * l.alpha[-1]
    l.alpha.adj[-1] <- l.alpha.adj[-1] + mm.u.diag * l.alpha[-T.tot]
    param.vec <- bin.counts + l.alpha.adj
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
    epsilon.out = epsilon.out, lambda.out = lambda.out
  ) 
}
