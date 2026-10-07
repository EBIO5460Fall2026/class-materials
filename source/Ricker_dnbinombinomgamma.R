# ---Negative binomial-binomial-gamma pmf for Ricker model
#
# Requires a standard number of increments and has no automatic resolution
# setting. Because of this, the whole vector of Ntp1 can be done at once and it
# is faster. See notes for Ricker_dnbinomgamma if higher accuracy is needed. The
# integral uses Eulers-McLaurin trapezoid geometry on P[R=r] - i.e. equal
# increments on the probability scale. This effectively gives us an adaptive
# step size - we take a smaller step size on the R scale where P[R=r] is high.
# This follows Numerical Recipes in C, trapzd.c and qtrap.c. The ptol and nmax
# values are quite course but sufficient for model fitting.
#
Ricker_dnbinombinomgamma <- function(Ntp1,Nt,p,R,alpha,kD,kE,log=FALSE,ptol=1e-3,nmax=6) {
# Integrate over R using an extended trapezoid rule
# Set initial increments
  Rl <- qgamma(ptol, shape = kE, scale = R / kE )
  Ru <- qgamma(1-ptol, shape = kE, scale = R / kE )
  Rhalf <- qgamma(0.5, shape = kE, scale = R / kE )
# One whole increment first and the halfway point
  psum <- 0.5 * ( Ricker_dnbinombinom_d(Ntp1,Nt,p,Rl,alpha,kD)
                  + Ricker_dnbinombinom_d(Ntp1,Nt,p,Ru,alpha,kD) )
  psum <- psum + Ricker_dnbinombinom_d(Ntp1,Nt,p,Rhalf,alpha,kD)
# Progressively smaller increments (halving size each time)
  for (n in 2:nmax) {
    totincs <- 2^n                #The total number of increments (old + new)
    newpts <- 2^(n-1)             #The number of new points
    nwdel <- (1-2*ptol)/newpts    #New delta
    p_evals <- seq(ptol+0.5*nwdel,1-ptol-0.5*nwdel,length.out=newpts)
    R_evals <- qgamma(p_evals, shape = kE, scale = R / kE )
    pmids <- 0                    #Initialize to zero
    for (r in R_evals) {
      pmids <- pmids + Ricker_dnbinombinom_d(Ntp1,Nt,p,r,alpha,kD)
    }
    psum <- psum + pmids          #Add the new points to previous ones
    pinc <- (1-2*ptol)/totincs    #Current increment size on probability scale
    pint <- psum * pinc           #Current approximation for p
  }
  if (log==FALSE) {
    return( pint )
  } # else
  return( log(pint) )
}

#----Negative binomial-binomial (demographic) pmf for Ricker model
#
# Translate to C? Can't avoid the loop, and this algorithm involves lots of
# memory allocation operations in R.
#
Ricker_dnbinombinom_d <- function(Ntp1,Nt,p,R,alpha,k,log=FALSE){
  sump <- rep(NA,length(Ntp1)) #Sum of the probabilities
  for (i in 1:length(Ntp1)){
    f <- 0:Nt[i] #f is number of females
    sump[i] <- sum( exp( dbinom(f,Nt[i],p,log=TRUE) +
        dnbinom(Ntp1[i],size=k*f+(f==0),
                mu=(1/p)*R*f*exp(-alpha*Nt[i]),log=TRUE) ) )
    #We use the log scale to ensure small probs multiply accurately
    #Add 1 to size parameter of NB when females=0 to avoid NaN's.
  }
  if (log==FALSE) {
    return( sump )
  } #else
  return( log(sump) )
}

#----Variance of the NB-binomial-gamma
#
Ricker_nbinombinomgamma.var <- function(Nt,R,alpha,kD,kE) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  dhvar <- Ricker(Nt,R,alpha)^2 / (kD*Nt)   #raw dhvar for no-sex model
  evar <- Ricker(Nt,R,alpha)^2 / kE         #env variance
  return( poisvar + sexvar + 2*dhvar + evar )
}

Ricker <- function(Nt, R, alpha) {
    mu = Nt * R * exp(-1 * alpha * Nt)
    return(mu)
}
