#' The Probability of Observing a Zero Value for a Tweedie Density
#' 
#' @description The probability that the variable takes the value of zero.
#'
#' @usage find_lambda(mu, phi, power)
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param power the power parameter \eqn{p} (sometimes denoted \eqn{\xi}{xi}).
#' 
#' @return The value of \eqn{\lambda}{lambda} when \eqn{1 < p < 2} such that \eqn{P(Y=0) = \exp(-\lambda)}{P(Y=0) = exp(-lambda)}. When \eqn{p>2}{power > 2}, a vector of zeros is returned.
#' 
#' @references 
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' @export
#' @examples
#' lambda <- find_lambda(mu = 1:3, phi = 1, power = 1.1)
#' exp( -lambda)
#' 
#' # When p > 2, there is zero probability that Y = 0:
#' lambda  <- find_lambda(mu = 1, phi = 1, power = 3.1)

find_lambda <- function(mu, phi, power){
  # Computes the value of lambda, such that P(Y = 0 ) = exp( -lambda) when 1 < p < 2
  if (power >= 2) {
    rep(0, length = length(mu) )
  } else {
    mu ^ (2 - power) / (phi * (2 - power) )
  }
}
