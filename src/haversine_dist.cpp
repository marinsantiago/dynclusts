// -----------------------------------------------------------------------------
// Matrix of Haversine distances
// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "haversine_dist.h"

using namespace Rcpp;


/**
 * Computes a matrix of Haversine distances based on a matrix of geographical  \
 * coordinates, where the first column corresponds to the longitude, while the \
 * second column corresponds to the latitude.                                  \
 *                                                                             \ 
 * @param coords Rcpp::NumericMatrix A matrix with two columns and n rows.     \
 *    The first column corresponds to the longitude, while the second column   \
 *    corresponds to the latitude.                                             \
 *                                                                             \ 
 * @return Rcpp::NumericMatrix A matrix of size n x n, with the Haversine      \
 * distances between the points in coords.                                     \
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix haver_dist(Rcpp::NumericMatrix & coords) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  int n = coords.nrow();
  Rcpp::NumericMatrix dist(n, n);
  const double R = 6378137.0; // Earth radius in meters
  
  // Convert degrees to radians ------------------------------------------------
  Rcpp::NumericVector lon(n);
  Rcpp::NumericVector lat(n);
  for (int i = 0; i < n; i++) {
    lon[i] = coords(i, 0) * M_PI / 180.0;
    lat[i] = coords(i, 1) * M_PI / 180.0;
  }
  
  // Compute distance ----------------------------------------------------------
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      double dlon = lon[j] - lon[i];
      double dlat = lat[j] - lat[i];
      double a = pow(sin(dlat / 2.0), 2.0) +
        cos(lat[i]) * cos(lat[j]) * pow(sin(dlon / 2.0), 2.0);
      double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
      double d = R * c;
      dist(i, j) = d;
      dist(j, i) = d; // symmetric
    }
  }
  
  return dist;
}
