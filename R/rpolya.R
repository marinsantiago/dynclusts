# ------------------------------------------------------------------------------
# Random number generation from the Polya(a, b) distribution 
# ------------------------------------------------------------------------------

# This function generates one random draw
rpolya <- function(a, b, trunc. = 200) {
  k <- 0:trunc.
  denom <- (a + k) * (b + k) / 2
  max(sum(rexp(trunc. + 1) / denom), 1e-10)
}
