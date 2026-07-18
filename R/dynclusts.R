#' Bayesian nonparametric modeling of dynamic clusters through an autoregressive 
#' logistic-beta Stirling-gamma (AR-LB-SG) process
#'
#' This function identifies latent dynamic clusters through an autoregressive 
#' logistic-beta Stirling-gamma (AR-LB-SG) process as described in 
#' Marin et al. (2026+).
#'
#' @param y A panel data matrix of size \eqn{T}-by-\eqn{n}, where each of the 
#'    \eqn{T} rows corresponds to a time point and each of the \eqn{n} columns 
#'    corresponds to a location.
#' @param coords A two column matrix with the geographical coordinates of each
#'    location. The first column corresponds to the longitude, while the second
#'    column corresponds to the latitude.
#' @param X A list of size \eqn{n}, containing matrices of size 
#'    \eqn{T}-by-\eqn{p} of covariates associated with each location. More
#'    precisely, each of the \eqn{T} rows contains \eqn{p} predictors at time 
#'    \eqn{t\in\{1,\dots,T\}} associated with the location 
#'    \eqn{i\in\{1,\dots,n\}}.
#' @param theta.0 Parameter \eqn{\theta_{0}} of the prior centering measure.
#' @param sigma2.0 Parameter \eqn{\sigma_{0}^{2}} of the centering measure.
#' @param a.0 Parameter \eqn{a_{0}} of the prior centering measure.
#' @param b.0 Parameter \eqn{b_{0}} of the prior centering measure.
#' @param a.alpha If \code{Sg = TRUE}, \code{a.alpha} corresponds to the \eqn{a}
#'    parameter in the Stirling-gamma prior on the concentration parameter of a
#'    Dirichlet process, as in Zito et al. (2026). If \code{Sg = FALSE}, 
#'    \code{a.alpha} corresponds to the \bold{shape} parameter in the gamma 
#'    prior on the concentration parameter of a Dirichlet process.
#' @param b.alpha If \code{Sg = TRUE}, \code{b.alpha} corresponds to the \eqn{b}
#'    parameter in the Stirling-gamma prior on the concentration parameter of a
#'    Dirichlet process, as in Zito et al. (2026). If \code{Sg = FALSE}, 
#'    \code{b.alpha} corresponds to the \bold{scale} parameter in the gamma 
#'    prior on the concentration parameter of a Dirichlet process.
#' @param a.phi \bold{Shape} parameter in the gamma prior on \eqn{\varphi}.
#'    Default is \code{0.1}.
#' @param b.phi \bold{Rate} parameter in the gamma prior on \eqn{\varphi}.
#'    Default is \code{0.1}.
#' @param a.rho \bold{Shape} parameter in the inverse-gamma prior on \eqn{\rho}.
#'    Default is \code{0.1}.
#' @param b.rho \bold{Scale} parameter in the inverse-gamma prior on \eqn{\rho}.
#'    Default is \code{0.1}.
#' @param a.tau \bold{Shape} parameter in the inverse-gamma prior on 
#'    \eqn{\tau^{2}}. Default is \code{0.1}.
#' @param b.tau \bold{Scale} parameter in the inverse-gamma prior on 
#'    \eqn{\tau^{2}}. Default is \code{0.1}.
#' @param sigma.phi_mh Initial magnitude of the proposed move in an adaptive 
#'    Metropolis-Hastings algorithm targeting the full conditional distribution
#'    of \eqn{\varphi}. Default is \code{sqrt(5.7)}.
#' @param sigma.psi_mh Initial magnitude of the proposed move in an adaptive 
#'    Metropolis-Hastings algorithm targeting the full conditional distribution
#'    of \eqn{\psi}. Default is \code{sqrt(0.001)}.   
#' @param H Truncation of the infinite mixture induced by a Dirichlet process 
#'    prior. Default is \code{25}.
#' @param Sg Logical. If \code{TRUE}, a Stirling-gamma prior is placed on the 
#'    concentration parameter of a Dirichlet process as in in 
#'    Zito et al. (20246). Otherwise, a gamma prior is employed. 
#'    Default is \code{TRUE}.
#' @param post.pred Logical. If \code{TRUE}, the posterior predictive 
#'    distribution of the model is stored and returned. Default is \code{TRUE}.
#' @param store.nu Logical. If \code{TRUE}, the posterior draws from the 
#'    parameter nu are stored and returned. Default is \code{TRUE}.    
#' @param store.epsilon Logical. If \code{TRUE}, the posterior draws from the 
#'    parameter epsilon are stored and returned. Default is \code{TRUE}.   
#' @param store.lambda Logical. If \code{TRUE}, the posterior draws from the 
#'    parameter lambda are stored and returned. Default is \code{TRUE}.   
#' @param verbose Logical. If \code{TRUE}, a progress bar is displayed on the
#'    console. Default is \code{TRUE}.
#' @param max.iters A positive integer corresponding to the total number of 
#'    MCMC iterations. Default is 20000.
#' @param burn.in A positive integer corresponding to the number of draws 
#'    discarded as burn-in. It should be smaller than \code{max.iters}.
#'    Default is 10000.
#' @param thin A positive integer specifying the period for saving samples.
#'    Default is 2.
#'
#' @return An object of S3 class, \code{"dynclust"}, containing:
#' \itemize{
#'   \item \code{theta.post}: An array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of 
#'   \eqn{\theta}.
#'   \item \code{sigma2.post}: An array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   \eqn{\sigma^{2}}.
#'   \item \code{omega.post}: An array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   the stick-breaking weights.
#'   \item \code{S.post}: An array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   the cluster allocations induced by the AR-LB-SG process.
#'   \item \code{K.post}: A matrix of size \code{n.draws}-by-\eqn{T} containing
#'   the posterior draws of the number of clusters at each time point.
#'   \item \code{alpha.post}: A matrix of size \code{n.draws}-by-\eqn{T}, 
#'   containing the posterior draws of \eqn{\boldsymbol{\alpha}}.
#'   \item \code{gamma.post}: A matrix of size \code{n.draws}-by-\eqn{n}, 
#'   containing the posterior draws of \eqn{\boldsymbol{\gamma}}.
#'   \item \code{psi.post}: A vector of size \code{n.draws}, containing the
#'   posterior draws of \eqn{\psi}.
#'   \item \code{varphi.post}: A vector of size \code{n.draws} containing the
#'   posterior draws of \eqn{\varphi}.
#'   \item \code{tau2.post}: A vector of size \code{n.draws} containing the
#'   posterior draws of \eqn{\tau^{2}}.
#'   \item \code{beta.post}: A matrix of size  \code{n.draws}-by-\eqn{p} 
#'   containing the posterior draws of \eqn{\boldsymbol{\beta}}.
#'   \item \code{rho2.post}: A vector of size \code{n.draws} containing the
#'   posterior draws of \eqn{\rho^{2}}.
#'   \item \code{nu.post}: If \code{store.nu = TRUE}, an array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   \eqn{\nu}.
#'   \item \code{epsilon.post}: If \code{store.epsilon = TRUE}, an array of size
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   \eqn{\epsilon}.
#'   \item \code{lambda.post}: If \code{store.lambda = TRUE}, an array of size
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the posterior draws of
#'   \eqn{\lambda}.
#'   \item \code{post.pred}: If \code{post.pred = TRUE}, an array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing draws from the posterior
#'   predictive distribution.
#'   \item \code{loglik}: An array of size 
#'   \code{T}-by-\code{n}-by-\code{n.draws} containing the log-likelihood of 
#'   each observation.
#'   \item Additional details of the function call.
#' }
#' 
#' @references
#'
#' S. Marin, B. Long,and A. H. Westveld (2026+), Bayesian nonparametric modeling
#' of dynamic pollution clusters through an autoregressive logistic-beta 
#' Stirling-gamma process. \emph{arXiv}, 2601.04625.
#' 
#' A. Zito, T. Rigon, and D. B. Dunson (20246), Bayesian Nonparametric Modeling 
#' of Latent Partitions via Stirling-Gamma Priors. \emph{Bayesian Analysis}, 21,
#' 139-166.
#' 
#' @author Santiago Marin
#' 
dynclust <- function(y, coords, X, theta.0, sigma2.0, a.0, b.0,
                     a.alpha, b.alpha, a.phi = 0.1, b.phi = 0.1, a.rho = 0.1,
                     b.rho = 0.1, a.tau = 0.1, b.tau = 0.1, 
                     sigma.phi_mh = sqrt(5.7), sigma.psi_mh = sqrt(0.001), 
                     H = 25L, Sg = TRUE, post.pred = TRUE, store.nu = TRUE, 
                     store.epsilon = TRUE, store.lambda = TRUE, verbose = TRUE,
                     max.iters = 20000L, burn.in  = 10000L, thin = 2L) {
  
  # Input validation -----------------------------------------------------------
  input.validation.dynclust(
    y = y, coords = coords, X = X, theta.0 = theta.0, sigma2.0 = sigma2.0,
    a.0 = a.0, b.0 = b.0, a.alpha = a.alpha, b.alpha = b.alpha, a.phi = a.phi,
    b.phi = b.phi, a.rho = a.rho, b.rho = b.rho, a.tau = a.tau, b.tau = b.tau,
    sigma.phi_mh = sigma.phi_mh, sigma.psi_mh = sigma.psi_mh, H = H, Sg = Sg,
    post.pred = post.pred, store.nu = store.nu, store.epsilon = store.epsilon,
    store.lambda = store.lambda, verbose = verbose, max.iters = max.iters,
    burn.in = burn.in, thin = thin
  )
  # Run MCMC algorithm ---------------------------------------------------------
  mcmc.out <- mcmc_sampler(
    y = y, coords = coords, X = X, theta.0 = theta.0, sigma2.0 = sigma2.0,
    a.0 = a.0, b.0 = b.0, a.alpha = a.alpha, b.alpha = b.alpha, a.phi = a.phi,
    b.phi = b.phi, a.rho = a.rho, b.rho = b.rho, a.tau = a.tau, b.tau = b.tau,
    sigma.phi_mh = sigma.phi_mh, sigma.psi_mh = sigma.psi_mh, H = H, Sg = Sg, 
    post.pred = post.pred, store.nu = store.nu, store.epsilon = store.epsilon,
    store.lambda = store.lambda, verbose = verbose, max.iters = max.iters,
    burn.in = burn.in, thin = thin
  )
  # Prepare the returns --------------------------------------------------------
  dims.y <- dim(y)
  mcmc.out[c("y", "X", "coords")] <- list(y, X, coords)
  mcmc.out[c("n.obs", "time.points")] <- list(dims.y[2], dims.y[1])
  mcmc.out["n.covariates"] <- ncol(X[[1]])
  mcmc.out[c("theta.0", "sigma2.0")] <- list(theta.0, sigma2.0)
  mcmc.out[c("a.0", "b.0", "a.rho", "b.rho")] <- list(a.0, b.0, a.rho, b.rho)
  mcmc.out[c("a.alpha", "b.alpha")] <- list(a.alpha, b.alpha)
  mcmc.out[c("a.tau", "b.tau")] <- list(a.tau, b.tau)
  mcmc.out[c("H", "Sg")] <- list(H, Sg)
  mcmc.out[c("max.iters", "burn.in", "thin")] <- list(max.iters, burn.in, thin)
  mcmc.out["n.draws"] <- length(mcmc.out$tau2)
  rm(dims.y); gc()
  class(mcmc.out) <- "dynclusts"
  mcmc.out
}
