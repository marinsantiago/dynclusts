#' Matrix of posterior co-clustering probabilities
#'
#' Computes the matrix of posterior co-clustering probabilities at a specific 
#' time point for an object of class \code{"dynclust"}.
#'
#' @param object An object of class \code{"dynclust"}.
#' @param time.point Integer indicating the time point at which we would like 
#'    to compute the co-clustering probabilities at. Default is \code{1}.
#' 
#' @return Matrix of posterior co-clustering probabilities
#'
#' @author Santiago Marin
#'
coclust_probs <- function(object, time.point = 1L) {
  if (!is.dynclust(object)) stop("object should be of class 'dynclusts'")
  # Pre-compute constants ------------------------------------------------------
  n <- object$n.obs
  dd <- object$n.draws
  T.tot <- object$time.points
  if (!(time.point %in% seq_len(T.tot))) stop("time.point out of range")
  # Compute co-clustering matrix -----------------------------------------------
  coclustering <- matrix(0, nrow = n, ncol = n)
  for (iter in seq_len(dd)) {
    # Clusters allocations at the corresponding time point
    cluster <- object$S.post[time.point,,iter]
    # Co-clusters
    coclustering <- coclustering + outer(cluster, cluster, FUN = "==")
  }
  coclustering.out <- coclustering / dd  # Normalize to get post. probabilities
  diag(coclustering.out) <- rep(1, n) # The diagonal is always one
  coclustering.out
}


#' Co-clustering matrix
#'
#' Computes the co-clustering matrix at a specific time point based on a matrix 
#' of cluster allocations. 
#'
#' @param S Matrix of cluster allocations.
#' @param time.point Integer indicating the time point at which we would like 
#'    to compute the co-clustering probabilities at. Default is \code{1}.
#' 
#' @return Co-clustering matrix
#'
#' @author Santiago Marin
#'
coclust_point <- function(S, time.point = 1L) {
  if (not.y(S)) stop("S must be a numeric matrix")
  if (!(time.point %in% seq_len(nrow(S)))) stop("time.point out of range")
  n <- ncol(S)
  # Compute co-clustering matrix -----------------------------------------------
  labels <- S[time.point,]
  coclusts <- outer(labels, labels, FUN = "==") * 1
  diag(coclusts) <- rep(1, n)  # The diagonal is always one
  coclusts
}


#' Matrix of lagged ARI values
#'
#' Computes the matrix of lagged ARI values based on a matrix of cluster 
#' allocations. 
#'
#' @param S Matrix of cluster allocations.
#' 
#' @return Matrix of lagged ARI values
#'
#' @author Santiago Marin
#'
lagged_ARI <- function(S) {
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Package 'salso' is required. Please install it.")
  }
  if (not.y(S)) stop("S must be a numeric matrix")
  T.tot <- nrow(S)
  # Compute lagged ARI matrix --------------------------------------------------
  lag.ari <- matrix(NA, nrow = T.tot, ncol = T.tot)
  for (t1 in seq_len(T.tot)) {
    for (t2 in seq_len(T.tot)) {
      lag.ari[t1, t2] <- salso::ARI(S[t1,], S[t2,])
    }
  }
  diag(lag.ari) <- rep(1, T.tot) # The diagonal is always one
  lag.ari
}
