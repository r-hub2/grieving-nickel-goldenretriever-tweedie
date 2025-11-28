rtweedie <- function(n, xi = NULL, mu, phi, power = NULL){
  
  ### BEGIN preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  out <- sort_notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!
  
  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  out <- check_inputs(y = n, 
                      mu = mu, 
                      phi = phi, 
                      power = power,
                      type = "random")
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )

  
  # IDENTIFY SPECIAL CASES
  out <- special_cases(n, mu, phi, power)
  f <- out$f
  special_p_cases <- out$special_p_cases
  special_y_cases <- out$special_y_cases

  ### END preliminary work
  
  
  if (power == 2) {
    alpha <- (2 - power) / (1 - power)
    gam <- phi * (power - 1) * mu ^ (power - 1)
    rt <- stats::rgamma( n, 
                  shape = 1 / phi, 
                  scale = gam )
  }
  
  if ( power > 2) {
    rt <- qtweedie( stats::runif(n),
                    mu = mu,
                    phi = phi, 
                    power = power)
  }
  
  if ( (power > 1) & (power < 2) ) {
    # Two options:  As above or directly.
    # Directly is faster
    rt <- array( dim = n, NA)
    
    lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
    alpha <- (2 - power) / (1 - power)
    gam <- phi * (power - 1) * mu ^ (power - 1)
    
    N <- stats::rpois(n, 
               lambda = lambda)
    for (i in (1:n) ){
      rt[i] <- stats::rgamma(1, 
                      shape = -N[i] * alpha, 
                      scale = gam[i])
    }
  }
  as.vector(rt)
}



