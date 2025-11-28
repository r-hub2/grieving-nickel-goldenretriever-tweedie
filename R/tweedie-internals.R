#' Tweedie internal functions
#'
#' @name Tweedie internals
#' @aliases dtweedie_logv_bigp dtweedie_logw_smallp dtweedie_interp dtweedie_jw_smallp dtweedie_kv_bigp dtweedie_Fortran dtweedie_Inversion_Report dtweedie_Inversion_Threemethods dtweedie_dldphi_saddle dtweedie_dlogfdphi dtweedie_logl dtweedie_dldphi dtweedie_logl_saddle dtweedie_series_bigp dtweedie_series_smallp logLiktweedie stored_grids twcomputation sort_notation check_inputs special_cases
#' @title Tweedie internal function
#' @description Internal tweedie functions. \bold{These are not to be called by the user.}
#'
#' @usage
#' dtweedie_dlogfdphi(y, mu, phi, power)
#' dtweedie_logl(phi, y, mu, power)
#' dtweedie_logl_saddle( phi, power, y, mu, eps=0)
#' dtweedie_logv_bigp( y, phi, power)
#' dtweedie_logw_smallp(y, phi, power)
#' dtweedie_interp(grid, nx, np, xix.lo, xix.hi,p.lo, p.hi, power, xix)
#' dtweedie_jw_smallp(y, phi, power )
#' dtweedie_kv_bigp(y, phi, power)
#' dtweedie_series_bigp(power, y, mu, phi)
#' dtweedie_series_smallp(power, y, mu, phi)
#' stored_grids(power)
#' check_inputs(y, mu, phi, power, type = "standard")
#' sort_notation(xi = NULL, power = NULL)
#' special_cases(y, mu, phi, power, type="PDF", verbose = FALSE)
#'
#'
#' @param y the vector of responses.
#' @param y the vector of responses.
#' @param power the value of \eqn{p}{power} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^p}{var(Y) = phi * mu^power}.
#' @param xi a synonym for \code{power}.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion \eqn{\phi}{phi}.
#' @param grid the interpolation grid necessary for the given value of \eqn{p}{power}.
#' @param nx the number of interpolation points in the \eqn{\xi}{xi} dimension.
#' @param np the number of interpolation points in the \eqn{p}{power}-dimension.
#' @param xix.lo the lower value of the transformed \eqn{\xi}{xi} value used in the interpolation grid. (Note that the value of \eqn{\xi}{xi} is from \eqn{0} to \eqn{\infty}{infty}, and is transformed such that it is on the range \eqn{0} to \eqn{1}.)
#' @param xix.hi the higher value of the transformed \eqn{\xi}{xi} value used in the interpolation grid.
#' @param p.lo the lower value of the \eqn{p} value used in the interpolation grid.
#' @param p.hi the higher value of the \eqn{p} value used in the interpolation grid.
#' @param xix the value of the transformed \eqn{\xi}{xi} at which a value is sought.
#' @param verbose logical; if \code{TRUE}, some details of the algorithm are returned.
#' @param type description
#' @param eps the offset in computing the variance function in the saddlepoint approximation. The default is \code{eps=1/6} (as suggested by Nelder and Pregibon, 1987).
#' @param type in \code{check_inputs}, the type of function for which inputs should be checked (one of \code{"standard"} (for \code{dtweedie} and \code{ptweedie}; the default), \code{"random"} (for \code{rtweedie}) or \code{"quantile"} (for \code{qtweedie})); in \code{special_cases}, one of \code{"PDF"} (the probability density function; the default) or \code{"CDF"} (the cumulative distribution function).
#'
#' @author Peter Dunn (\email{pdunn2@usc.edu.au})
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008).
#' Evaluation of Tweedie exponential dispersion model densities by Fourier inversion.
#' \emph{Statistics and Computing}, 
#' \bold{18}, 73--86.
#' \doi{10.1007/s11222-007-9039-6}
#' 
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 
#' Jorgensen, B. (1997).
#' \emph{Theory of Dispersion Models}.
#' Chapman and Hall, London.
#'
#' 	Nelder, J. A. and Pregibon, D. (1987).
#' 	An extended quasi-likelihood function
#' 	\emph{Biometrika},
#' 	\bold{74}(2), 221--232.
#' 	\doi{10.1093/biomet/74.2.221}
#' 	
#' @keywords internal


dtweedie_dldphi_saddle <- function(phi, mu, power, y){
  # Calculates the derivative of log f wrt phi
  # where the density is the saddlepoint density
  
  # Peter Dunn
  # 13 August 2002
  
  dev <- tweedie_dev( power = power, 
                      y = y, 
                      mu = mu)
  l <-  (-1) / (2 * phi) + dev / (2 * phi ^ 2)
  
  -2 * sum(l)
}




#############################################################################


dtweedie_logl <- function(phi, y, mu, power) {
  # Computes the log-likelihood for
  # a Tweedie density.  
  
  # Peter Dunn
  # 26 April 2001
  
  sum( log( dtweedie( y = y, 
                      mu = mu, 
                      phi = phi, 
                      power = power) ) )
  
}


#############################################################################


logLiktweedie <- function(glm.obj, dispersion = NULL) {
  # Computes the log-likelihood for
  # a Tweedie glm.  
  
  # Peter Dunn
  # 19 October 2017
  
  p <- get("p", envir = environment(glm.obj$family$variance))
  if (p == 1) message("*** Tweedie index power = 1: Consider using  dispersion = 1  in call to  logLiktweedie().\n")
  
  AICtweedie(glm.obj, 
             dispersion = dispersion, 
             k = 0, 
             verbose = FALSE) / (-2)
}


#############################################################################

dtweedie_logl_saddle <- function( phi, power, y, mu, eps=0){
  # Calculates the log likelihood of Tweedie densities
  # where the density is the saddlepoint density
  
  # Peter Dunn
  # 01 May 2001
  sum( log( dtweedie_saddle(power = power, 
                            phi = phi, 
                            y = y, 
                            mu = mu,
                            eps = eps) ) )
  
}




#############################################################################


dtweedie_series_bigp <- function(power, y, mu, phi){ 
  
  # 
  # Peter K Dunn 
  # 02 Feb 2000 
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2) stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  result <- dtweedie_logv_bigp(power = power, 
                               y = y, 
                               phi = phi)
  logv <- result$logv
  
  theta <- mu ^ (1 - power) / ( 1 - power )
  kappa <- mu ^ (2 - power) / ( 2 - power )
  
  logfnew <- (y * theta - kappa) / phi - log( pi * y) + logv
  f <- exp( logfnew )
  
  list(density = f, 
       logv = logv, 
       lo = result$lo, 
       hi = result$hi )
  
}


#############################################################################

dtweedie_dldphi <- function(phi, mu, power, y ){
  # Calculates the log-likelihood
  # function, wrt phi, for p>2.  In particular, it returns
  #    sum{ d(log f)/d(phi) } = d( log-likelihood )/d(phi).
  # The mle of phi can be found, therefore, by setting this to zero.
  # y  is generally a vector of observed values.
  
  #
  # Peter Dunn
  # 31 Jan 2001
  #
  
  if ( (power != 2 ) & ( power != 1 ) ) {
    
    k <- phi ^ (1 / (power - 2))
    # cat("k=",k,"\n")
    # cat("phi=",phi,"\n")
    # cat("power=",power,"\n")
    if ( k < 1  & k > 0 ) {
      # Use the transform f(y; mu, phi) = c f(c*y; c*mu, c^(2-p)*phi)
      # and differentiate with c=phi^(1/(p-2)):
      #    d log f / d phi = c^(2-p) * {df(cy; c*mu, 1)/dphi} / f(cy; c*mu, 1)
      f <- dtweedie( y = k * y, 
                     power = power, 
                     mu = k * mu, 
                     phi = 1 )
      d <- dtweedie_dlogfdphi( y = k * y, 
                               power = power, 
                               mu = k * mu, 
                               phi = 1 )
      # Note:  We need dlogf/dphi = dlogf.dphi * f
      top <- d * f
      d <- -2* sum( top / f * k ^ (2 - power) )
      
    } else{
      # Compute directly
      d <- -2 * sum( dtweedie_dlogfdphi(y = y, 
                                        power = power, 
                                        mu = mu, 
                                        phi = phi) )
    }
  } else{
    # Cases p == 1 and  p == 2 
    d <- -2 * sum( dtweedie_dlogfdphi(y = y, 
                                      power = power, 
                                      mu = mu, 
                                      phi = phi) )
  }
  d
}



#############################################################################

dtweedie_dlogfdphi <- function(y, mu, phi, power)
{
  #
  # Calculates d(log f)/d(phi) for the Tweedie
  # densities.
  # It is used, for example, in mle fitting of phi.  We would then
  # sum over  y  and set this function to 0.
  #
  #
  # Peter Dunn
  # 31 Jan 2001
  #
  
  p <- power
  a <- (2 - p) / (1 - p)
  if(length(phi) == 1) {
    phi <- array(dim = length(y), phi)
  }
  if(length(mu) == 1) {
    mu <- array(dim = length(y), mu)
  }
  A <- (y * mu ^ (1 - p)) / (phi ^ 2 * (p - 1))
  B <- mu^(2 - p)/(phi ^ 2 * (2 - p))
  
  if(power > 2) {
    f <- array(dim = c(length(y)))
    # Here, we evaluate logv and everything as normal.
    # If logv has infinite values then we resort to other tactics.
    
    kv <- dtweedie_kv_bigp(power = power, 
                           phi = phi, 
                           y = y)$kv
    dv.dphi <- (kv * (a - 1)) / phi
    out.logv <- dtweedie_logv_bigp(power = power, 
                                   phi = phi, 
                                   y = y)
    
    # Now see if this causes problems.
    logv <- out.logv$logv
    
    # Now detect problem computing  logv  and remedy them
    probs <- (is.infinite(logv)) | (is.nan(logv)) | (y < 1)
    
    if(any(probs)) {
      
      # OK then:  Troubles computing log(V).
      # best we can do is use definition I think.
      delta <- 1.0e-5
      a1 <- dtweedie(power = power, 
                     phi = phi[probs], 
                     mu = mu[probs], 
                     y = y[probs])
      a2 <- dtweedie(power = power, 
                     phi = phi[probs] + delta, 
                     mu = mu[probs], 
                     y = y[probs])
      f[probs] <- (log(a2) - log(a1) ) / delta
      
    }
    
    f[!probs] <- A[!probs] + B[!probs] + dv.dphi[ !probs] / exp(logv[!probs])
  }
  #  END p>2
  
  if(power == 2) {
    f <-   -log(y)  + ( y / mu ) + digamma(1 / phi) - 1 + log( mu * phi )
    f <- f / (phi ^ 2)
  }
  
  if(power == 1) {
    
    f <- mu - y - y * log(mu / phi) + y * digamma(1 + (y / phi))
    f <- f / (phi ^ 2)
  }
  
  if((power > 1) && (power < 2)) {
    # We need to treat y==0 separately.
    # We fill  f  with the Y=0 case
    f <- array(  dim = length(y), 
                 mu ^ (2 - power) / ( phi ^ 2 * (2 - power) )  )
    jw <- dtweedie_jw_smallp(power = power, 
                             phi = phi[y > 0], 
                             y = y[y > 0])$jw
    dw.dphi <- (jw * (a - 1)) / phi[y > 0]
    logw <- dtweedie_logw_smallp(power = power, 
                                 phi = phi[y > 0], 
                                 y = y[y > 0])$logw
    f[y>0] <- A[y > 0] + B[y > 0] + dw.dphi / exp(logw)
  }
  f
}




