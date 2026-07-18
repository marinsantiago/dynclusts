#' Clustering estimates
#'
#' Identifies the partition that minimizes an user-specified distance function 
#' for an object of class \code{"dynclust"}.
#'
#' @param object An object of class \code{"dynclust"}.
#' @param dist A character string denoting which distance function should be 
#'    used. Options are either \code{"VI"} for the  variation of information 
#'    oss function or \code{"Binder"} for Binder's loss function. Default is 
#'    \code{"VI"}
#' @param return.clusters.means Logical. Whether the means of each cluster 
#'    should be returned. Default is \code{TRUE}. 
#'
#' @return Clustering of the observations.
#'
#' @author Santiago Marin
#'
clustering <- function(object, dist = "VI", return.clusters.means = TRUE) {
  if (!requireNamespace("BNPmix", quietly = TRUE)) {
    stop("Package 'BNPmix' is required. Please install it.")
  }
  # Input validation -----------------------------------------------------------
  if (!is.dynclust(object)) stop("object should be of class 'dynclusts'")
  check <- not.logic(return.clusters.means)
  if (check) stop("'return.clusters.means' must be logical") 
  check <- !(dist %in% c("VI", "Binder"))
  if (check) stop("'dist' must be either 'VI' or 'Binder'") 
  # Extract constants ----------------------------------------------------------
  n <- object$n.obs
  dd <- object$n.draws
  T.tot <- object$time.points
  # Clustering at each time point  ---------------------------------------------
  clusters.out <- matrix(nrow = T.tot, ncol = n)
  if (return.clusters.means) theta.out <- matrix(nrow = T.tot, ncol = n)
  for (tt in seq_len(T.tot)) {
    clusters.draws <- matrix(nrow = dd, ncol = n)
    for (iter in seq_len(dd)) clusters.draws[iter, ] <- object$S.post[tt, ,iter]
    #bnpmix.object <- list(clust = clusters.draws - 1)
    bnpmix.object <- list(clust = clusters.draws)
    class(bnpmix.object) <- "BNPdens"
    BNPmix.out <- BNPmix::partition(bnpmix.object, dist = dist)
    min.lower.bound.expected.loss <- which.min(BNPmix.out$scores)
    clusters.t <- BNPmix.out$partitions[min.lower.bound.expected.loss,]
    clusters.out[tt,] <- clusters.t
    if (return.clusters.means) {
      clusters.idx <- unique(clusters.t)
      theta.t <- matrix(nrow = dd, ncol = n)
      for (iter in seq_len(dd)) theta.t[iter, ] <- object$theta.post[tt, ,iter]
      # Iterate over clusters
      for (k in seq_along(clusters.idx)) {
        current.obs.clust.k <- clusters.t == clusters.idx[k]
        theta.clust.k <- mean(theta.t[,current.obs.clust.k])
        theta.out[tt, current.obs.clust.k] <- theta.clust.k
      }
    }
  }
  out <- list(clusters = clusters.out)
  if (return.clusters.means) out[["clusters.means"]] <- theta.out
  out
}
