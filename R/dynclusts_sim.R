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
#' @param pi.stay Probability that a given station will remain in the same 
#'    cluster from one time point to another. Default is \code{0.3}.
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
                          varphi = 100, sigma2 = 1, pi.stay = 0.3, 
                          balanced = FALSE,spatial.effect = TRUE, 
                          covariates.effect = TRUE) {
  
  # Input validation -----------------------------------------------------------
  input.validation.dynclust.sim(
    n.obs, time.points, n.clusts, n.covariates, country, tau2, rho2, varphi, 
    sigma2, pi.stay, balanced, spatial.effect, covariates.effect
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
  # Cluster probabilities
  if (balanced) {
    # All clusters have equal probability
    p.equal <- rep(1 / n.clusts, n.clusts)
  } else {
    # Probability of the largest cluster
    p.large <- 0.7 
    # Probability of the remaining clusters
    p.rest <- rep((1 - 0.7) / (n.clusts - 1), n.clusts - 1)
  }
  # Clusters at time t = 1
  probs <- if (balanced) p.equal else c(p.rest, p.large)
  S[1,] <- int_sampling(n.clusts, n.obs, probs)
  theta[1,] <- theta.candidates[S[1,]]
  # Generate y
  y[1,] <- rnorm(n.obs, linpred[1,] + theta[1,], sigma2[1,])
  # Clusters at times t = 2, 3, ..., T
  for (tt in 2:time.points) {
    # Whether the data points should stay in their current clusters
    clust.stay <- rbinom(n.obs, 1, pi.stay)
    # Iterate over data points
    for (ii in seq_len(n.obs)) {
      if (clust.stay[ii] == 1) {
        S[tt, ii] <- S[tt - 1, ii] # Remain in the same cluster
      } else {
        # Update the probabilities of each cluster
        probs <- if (!balanced) {
          if (tt %% 2 == 0) c(p.large, p.rest) else c(p.rest, p.large)
        } else { p.equal }
        # Jump to a new cluster with the updated cluster probabilities
        S[tt, ii] <- int_sampling(n.clusts, 1, probs)
      }
    }
    theta[t,] <- theta.candidates[S[t,]]
    # Generate y
    y[t,] <- rnorm(n.obs, linpred[t,] + theta[t,], sigma2[t,])
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
