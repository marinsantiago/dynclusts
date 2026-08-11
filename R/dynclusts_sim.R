#' Simulate dynamic clusters 
#'
#' This function simulates dynamic clusters as in Marin et al. (2026+).
#'
#' @param n.obs Number of locations
#' @param time.points Number of time points.
#' @param n.clusts Number of clusters. Default is \code{3}.
#' @param n.covariates Number of covariates. Default is \code{5}.
#' @param country Country to pick locations from. Default is "United Kingdom".
#' @param tau2 Variance of the location-specific random effects. Default is 
#'    \code{2}.
#' @param rho2 Variance of the vector of coefficients, \eqn{\boldsymbol{\beta}}.
#'    Default is \code{5}.
#' @param varphi Strength of the spatial correlation between two locations.
#'    Default is \code{100}.
#' @param sigma2 Sampling variance, assumed to be common across all clusters. 
#'    Default is \code{1}.
#' @param prob.jump Probability of jumping, \eqn{\pi_{\text{jump}}}. Default is 
#'    \code{0.2}.
#' @param prob.pow Power to which \code{prob.jump} is elevated to. Default is
#'    \code{0.25}.
#' @param balanced logical. Whether cluster sizes should be balanced or not. 
#' @param spatial.effect logical. Whether spatial effects should be
#'    incorporated.
#' @param covariates.effect logical. Whether covariate effects should be 
#'    included.
#'    
#' @return An object of S3 class, \code{"dynclust.sim"}, containing:
#' \itemize{
#'   \item \code{y}: A panel data matrix of size \eqn{T}-by-\eqn{n}, where each 
#'   of the \eqn{T} rows corresponds to a time point and each of the \eqn{n} 
#'   columns corresponds to location.
#'   \item \code{S}: A matrix of size \eqn{T}-by-\eqn{n} with the cluster
#'   allocations.
#'   \item \code{coords}: A two column matrix with the geographical coordinates 
#'   of each location. The first column corresponds to the longitude, while the
#'   second column corresponds to the latitude.
#'   \item \code{theta}: A matrix of size \eqn{T}-by-\eqn{n} with each of the
#'   \eqn{\theta_{ti}} parameters.
#'   \item \code{sigma2}: A matrix of size \eqn{T}-by-\eqn{n} with each of the
#'   \eqn{\sigma_{ti}^{2}} parameters.
#'   \item \code{gamma}: If \code{spatial.effect = TRUE}, a vector of size 
#'   \eqn{n} with the location-specific random effects, \eqn{\gamma_{i}}.
#'   \item \code{K}: If \code{spatial.effect = TRUE}, the squared exponential
#'   kernel.
#'   \item \code{beta}: If \code{covariates.effect = TRUE}, a vector of size
#'   \code{n.covariates} with the coefficients associated with each of 
#'   covariate.
#'   \item \code{X}: If \code{covariates.effect = TRUE}, a list of size \eqn{n},
#'   containing matrices of size \eqn{T}-by-\eqn{p} of covariates associated
#'   with each location. More precisely, each of the \eqn{T} rows contains 
#'   \eqn{p} predictors at time \eqn{t\in\{1,\dots,T\}} associated with the 
#'   location \eqn{i\in\{1,\dots,n\}}.
#' }
#' 
#' @references
#'
#' S. Marin, B. Long,and A. H. Westveld (2026+), Bayesian nonparametric modeling
#' of dynamic pollution clusters through an autoregressive logistic-beta 
#' Stirling-gamma process. \emph{arXiv}, 2601.04625.
#' 
#' @author Santiago Marin
#'
dynclusts_sim <- function(n.obs, time.points, n.clusts = 3L, n.covariates = 5L,
                          country = "United Kingdom", tau2 = 2, rho2 = 1, 
                          varphi = 100, sigma2 = 1, prob.jump = 0.2, 
                          prob.pow = 0.25,  balanced = TRUE, 
                          spatial.effect = TRUE, covariates.effect = TRUE) {
  
  # Input validation -----------------------------------------------------------
  input.validation.dynclust.sim(
    n.obs = n.obs, time.points = time.points, n.clusts = n.clusts,
    n.covariates = n.covariates, tau2 = tau2, rho2 = rho2, varphi = varphi,
    sigma2 = sigma2, prob.jump = prob.jump, prob.pow = prob.pow, 
    balanced = balanced, spatial.effect = spatial.effect, 
    covariates.effect = covariates.effect
  )
  # Set the country and sample geo-points at random (on-land!) -----------------
  region_sf <- ne_countries_medium[ne_countries_medium$name == country, ]
  # NB: "ne_countries_medium" is stored in "sysdata.rda"
  # Bounding box
  bbox <- sf::st_bbox(region_sf)
  # "Special" countries with insular territories far away from the mainland. 
  # Consider only "continental" land. Others can be added later if needed.
  special.bounds <- if (country == "Chile") {
    c(xmin = -75.6493, ymin = -56.5375, xmax = -66.9754, ymax = -17.4983)
  } else if (country == "United States of America") { # 48-cont. states
    c(xmin = -124.7331, ymin = 24.4467, xmax = -66.950, ymax = 49.3845)
  } else if (country == "New Zealand") {
    c(xmin = 166.2199, ymin = -47.3178, xmax = 178.5362, ymax = -34.2635)
  } else if (country == "France") {
    c(xmin = 166.2199, ymin = -47.3178, xmax = 178.5362, ymax = -34.2635)
  } else { NA }
  if (any(!is.na(special.bounds))) bbox <- sf::st_bbox(special.bounds)
  points <- list()
  count <- 0
  while (length(points) < n.obs) {
    lon <- runif(1, bbox["xmin"], bbox["xmax"])
    lat <- runif(1, bbox["ymin"], bbox["ymax"])
    pt <- sf::st_point(c(lon, lat)) |> sf::st_sfc(crs = sf::st_crs(region_sf))
    if (sf::st_intersects(pt, region_sf, sparse = FALSE)[1]) {
      points[[length(points) + 1]] <- pt
    }
  }
  geopoints.out <- sf::st_as_sf(do.call(c, points))
  # Matrix of geographical coordinates
  coords <- do.call(rbind, lapply(points, `[[`, 1L))
  colnames(coords) <- c("longitude", "latitude")
  # Generate covariates --------------------------------------------------------
  X <- replicate(
    n.obs, matrix(
      if (covariates.effect) runif(time.points * n.covariates) else 0,
      nrow = time.points, ncol = n.covariates
    ), simplify = FALSE
  )
  # Generate location-specific random effects ----------------------------------
  gamma. <- if (spatial.effect) {
    D <- (haver_dist(coords) / 1000)^2 # Squared distances (in Km)
    K <- make_posdef(tau2 * exp(-D / (2 * varphi^2)))  # Squared exp. kernel
    as.vector(rnorm(n.obs) %*% chol(K)) + rep(3, n.obs)
  } else { rep(0, n.obs) }
  # Generate coefficient vector and linear predictor ---------------------------
  beta. <- if (covariates.effect) {
    rnorm(n.covariates, mean = 3, sd = sqrt(rho2))
  } else { rep(0, n.covariates) }
  linpred <- do.call(rbind, X) %*% beta.
  dim(linpred) <- c(time.points, n.obs)
  linpred <- linpred + rep(gamma., each = time.points) # Add theta later
  # Generate theta, sigma2, and y ----------------------------------------------
  theta.candidates <- round(seq(6, 60, length.out = n.clusts))
  theta <- y <- matrix(nrow = time.points, ncol = n.obs)
  # Assumed same sigma2 across clusters
  sigma2 <- matrix(sigma2, nrow = time.points, ncol = n.obs)
  S <- matrix(nrow = time.points, ncol = n.obs) # Cluster allocations
  clusters <- seq_len(n.clusts)
  if (balanced) {
    # Clusters at time t = 1
    probs <- rep(n.obs / n.clusts, n.clusts)
    probs <- probs / sum(probs)
    S[1,] <- int_sampling(n.clusts, n.obs, probs)
    theta[1,] <- theta.candidates[S[1,]]
    # Generate y
    y[1,] <- rnorm(n.obs, linpred[1,] + theta[1,], sigma2[1,])
    # Times t = 2, 3, ..., T
    for (t in 2:time.points) {
      for (k in clusters) {
        idx.k <- which(S[t-1,] == k)
        n.k <- length(idx.k)
        S[t, idx.k] <- sample(
          c(k, setdiff(clusters, k)), size = n.k, replace = TRUE,
          prob = c(1 - prob.jump, prob.jump / 2, prob.jump / 2)
        )
      }
      theta[t,] <- theta.candidates[S[t,]]
      # Generate y
      y[t,] <- rnorm(n.obs, linpred[t,] + theta[t,], sigma2[t,])
    }
  } else {
    # Clusters at time t = 1
    big.clust.idxs <- sort(sample.int(n.obs, round(n.obs * 0.7), replace = F))
    # Choose at random which cluster is going to be the biggest one
    big.clust.k <- sample.int(n.clusts, 1L)
    S[1, big.clust.idxs] <- big.clust.k
    # Remaining clusters
    other.clusts.idxs <- setdiff(seq_len(n.obs), big.clust.idxs)
    other.clusts.k <- setdiff(clusters, big.clust.k)
    other.clusts.n <- length(other.clusts.idxs)
    S[1, other.clusts.idxs] <- sample(
      other.clusts.k, size = other.clusts.n, replace = TRUE
    )
    theta[1,] <- theta.candidates[S[1,]]
    # Generate y
    y[1,] <- rnorm(n.obs, linpred[1,] + theta[1,], sigma2[1,])
    # Times t = 2, 3, ..., T
    for (t in 2:time.points) {
      # Identify the biggest cluster at time t-1
      big.clust.k <- which.max(tabulate(S[t-1,]))
      big.clust.idxs <- which(S[t-1,] == big.clust.k)
      big.clust.n <- length(big.clust.idxs)
      # Randomly select the new biggest cluster 
      other.clusts.k <- setdiff(clusters, big.clust.k)
      big.clust.k.new <- sample(other.clusts.k, 1L)
      # Randomly select which observations are going to jump from the
      # current biggest cluster to the new one
      big.jump <- as.logical(rbinom(big.clust.n, 1, prob.jump ** prob.pow))
      big.jump.idx <- big.clust.idxs[big.jump]
      big.remain.idx <- big.clust.idxs[!big.jump]
      S[t, big.jump.idx] <- big.clust.k.new
      S[t, big.remain.idx] <- big.clust.k
      # Allocate the remaining observations
      for (k in other.clusts.k) {
        idx.k <- which(S[t-1,] == k)
        n.k <- length(idx.k)
        S[t, idx.k] <- sample(
          c(k, setdiff(clusters, k)), size = n.k, replace = TRUE,
          prob = c(1 - prob.jump, prob.jump / 2, prob.jump / 2)
        )
      }
      theta[t,] <- theta.candidates[S[t,]]
      # Generate y
      y[t,] <- rnorm(n.obs, linpred[t,] + theta[t,], sigma2[t,])
    }
  }
  # Returns --------------------------------------------------------------------
  out <- list(
    y = y, S = S, coords = coords, geopoints.out = geopoints.out, n.obs = n.obs,
    country = country, theta = theta, sigma2 = sigma2, time.points = time.points
  )
  if (spatial.effect) {
    out[["gamma"]] = gamma.
    out[["K"]] = K
    out[["varphi"]] = varphi
    out[["tau2"]] = tau2
  }
  if (covariates.effect) {
    out[["beta"]] = beta.
    out[["rho2"]] = rho2
    out[["X"]] = X
    out[["n.covariates"]] = n.covariates
  }
  class(out) <- "dynclust.sim"
  out
}
