#'  Series Evaluation for the Tweedie Distribution Function
#'
#' @description
#' Evaluates the distribution function (\acronym{df}) for Tweedie distributions 
#' with \eqn{1 < p < 2}{1 < p < 2}
#' using an infinite series, for given values of the dependent variable \code{y}, 
#' the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be in the case of evaluation problems.
#'
#' @usage ptweedie_series(q, power, mu, phi, verbose = FALSE, details = FALSE)
#' 
#' @param q vector of quantiles.
#' @param power the power parameter \eqn{p}{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param verbose logical; if \code{TRUE}, displays some internal computation details. The default is \code{FALSE}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution function and some details.
#' 
#' @return A numeric vector of densities.
#' 
#' @references
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' @examples
#' # Plot a Tweedie distribution function
#' y <- seq(0, 5, length = 100)
#' Fy <- ptweedie_series(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, Fy, type = "l", lwd = 2, ylab = "Distribution function")
#' 
#' @importFrom stats dpois 
#'
#' @keywords distribution
#'
#' @export
#' 
#' @aliases ptweedie.series

ptweedie_series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE) {
  ### NOTE: No notation checks

  # SET UP
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  drop <- 39

  # FIND THE LIMITS ON N, the summation index
  # The *lower* limit on N
  lambda <- max(lambda )
  logfmax <-  -log(lambda)/2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( ( estlogf > (logfmax - drop) ) & ( N > 1 ) ) {
    N <- max(1, N - 2)
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  lo.N <- max(1, floor(N) )
  
  
  # The *upper* limit on N
  lambda <- min( lambda )
  logfmax <-  -log(lambda) / 2
  estlogf <- logfmax
  N <- max( lambda )
  
  while ( estlogf > (logfmax - drop) ) {
    N <- N + 1
    estlogf <- -lambda + N * ( log(lambda) - log(N) + 1 ) - log(N)/2
  }
  hi.N <- max( ceiling(N) )
  if (verbose) cat("Summing over", lo.N, "to", hi.N, "\n")
  
  # EVALUATE between limits of N
  cdf <- array( dim = length(q), 0 )
  
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  tau    <- phi * (power - 1) * mu ^ ( power - 1 )
  alpha  <- (2 - power) / (1 - power)
  
  
  for (N in (lo.N : hi.N)) {
    # Poisson density
    pois.den <- dpois( N, lambda)
    
    # Incomplete gamma
    incgamma.den <- stats::pchisq(2 * q / tau, 
                           -2 * alpha * N )
    
    # What we want
    cdf <- cdf + pois.den * incgamma.den
    
  }
  
  cdf <- cdf + exp( -lambda )
  its <- hi.N - lo.N + 1
  
  if (details) {
    return( list( cdf = cdf,
                  iterations = its) )
  } else {
    return(cdf)
  }
  
}

#' @export
ptweedie.series <- function(q, power, mu, phi, verbose = FALSE, details = FALSE){ 
  .Deprecated("ptweedie_series", package = "tweedie")
  ptweedie_series(q, power, mu, phi, verbose = FALSE, details = FALSE)
}

