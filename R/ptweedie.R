ptweedie <- function(q, xi = NULL, mu, phi, power = NULL, verbose = FALSE) {

  ### BEGIN preliminary work

  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out <- sort_notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  
  if (verbose) cat("- Checking, resizing inputs\n")
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_inputs(q, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )

  # IDENTIFY SPECIAL CASES
  special_y_cases <- rep(FALSE, length(q))
  if (verbose) cat("- Checking for special cases\n")
  out <- special_cases(q, mu, phi, power, 
                       type = "CDF")
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases

  if (verbose & special_p_cases) cat("  - Special case for p used\n")
  if ( any(special_y_cases) ) {
    special_y_cases <- out$special_y_cases  
    if (verbose) cat("  - Special cases for first input found\n")
    f <- out$f
  }

  ### END preliminary work

  if ( special_p_cases ) {
    f <- out$f
  } else {
    # NOT special p case; ONLY special y cases 
    
    if ( power > 2 ) {
      # For p > 2 the only option is the inversion
      if ( any(!special_y_cases)) { 
        if (verbose) cat("- With p > 2: use inversion\n")
  
        f_TMP <- ptweedie_inversion(power   = power,
                                    q       = q[!special_y_cases],
                                    mu      = mu[!special_y_cases],
                                    phi     = phi[!special_y_cases],
                                    verbose = FALSE,
                                    details = FALSE)
        f[!special_y_cases] <- f_TMP
      }
    } else {
      # CASE 1 < p < 2
      # For 1 < p < 2, the two options are the series or inversion.
      # We avoid the series when p is near one, otherwise it is fine.
      
      # A few caveats.  Gustaov noted this case:
      # ptweedie(q=7.709933e-308, mu=1.017691e+01, phi=4.550000e+00, power=1.980000e+00)
      # which fails for the inversion, but seems to go fine for the series.
      # So the criterion should be a bit more detailed that just p<1.7...
      # But what?
      # Shallow second derivative of integrand is OK in principle... but hard to ascertain.
      
      # In a bug report by Gustavo Lacerda (April 2017), 
      # it seems that the original code here was a bit too 
      # harsh on the series, and a bit forgiving on the inversion.
      # Changed April 2017 to use the series more often.
      
      # ### OLD CODE:
      # if ( (power>1) & (power<2) ) {
      #     if ( power <1.7 ) {
      #        f <- ptweedie_series(power=power, q=y, mu=mu, phi=phi )
      #     } else{
      #        f <- ptweedie_inversion( power=power, q=y, mu=mu, phi=phi)
      #     }
      # }
    
      ### REVISED CODE:
      ### The choice of  1.999 is arbitrary.  Probably needs serious attention to decide properly
      ### Changed early 2017; thanks to Gustavo Lacerda
      #if ( power < 1.999) { 
      #  #### XXXXXXXXXXXXXXXXXXXXXXXXX This is arbitrary, and needs a closer look
      #  if (verbose) cat("- With 1 < p < 2: use series")
      #  
      #  f_TMP <- ptweedie_series(power = power, 
      #                           q     = y[!special_y_cases], 
      #                           mu    = mu[!special_y_cases], 
      #                           phi   = phi[!special_y_cases] )
      #  if (details) {
      #    f[!special_y_cases] <- f_TMP$cdf
      #    regions[!special_y_cases] <- f_TMP$regions
      #  } else {
      #    f[!special_y_cases] <- f_TMP
      #  }
      #} else{
        
      if ( any(!special_y_cases)) {
        if (verbose) cat("- With 1 < p < 2: use inversion TEMPORARILY")
        f_TMP <- ptweedie_inversion(power   = power,
                                    q       = q[!special_y_cases], 
                                    mu      = mu[!special_y_cases], 
                                    phi     = phi[!special_y_cases],
                                    verbose = FALSE,
                                    details = FALSE)
        f[!special_y_cases] <- f_TMP
      }
    }  
  }

  # Sanity fixes
  f <- as.vector(f)
  f[ f < 0 ] <- rep(0, sum(f < 0) )
  f[ f > 1 ] <- rep(1, sum(f > 1) )

  return(f)
}

