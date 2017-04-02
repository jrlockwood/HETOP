##kat updated to return cut estimates as well 
est_deref <- function(mprime, sprime, cutests, ngs, geticc, biascor = FALSE) {
  #   function to "standardize" the HETOP estimates as in our paper
  #   mus=vector of mean estimates
  #   sigs=vector of sigma estimates
  #   ngs=vector of group sample sizes
  #   cutests=estimated cuts
  #   biascor: TRUE=formulas in paper; FALSE=simpler standardizing as with jags
  #   will need to think about which identification is used
  #   for now, should not matter too much because we have equally sized groups
  
  N <- sum(ngs)
  G <- length(ngs)
  ntilde <- 1/sum((1/(ngs-1))*(1/G))
  p <- ngs/N
  
  if (biascor) {
    ssw_pw <- sum( (ntilde/(1+ntilde)) * p * sprime^2 )
    ssb_pw <- sum( p* (mprime - sum(p * mprime) )^2 ) - sum( (ntilde/(1+ntilde)) * (1/N) * (1-p) * sprime^2 )
    
    icchat_pw <- ssb_pw / (ssb_pw + ssw_pw)
    
    mstar <- (mprime - sum(p*mprime)) / sqrt(ssb_pw + ssw_pw)
    sstar <- sprime / sqrt(ssb_pw + ssw_pw)
    cutstar <- (cutests-sum(p*mprime))/sqrt(ssb_pw + ssw_pw)
    
    if (geticc == TRUE) {
      return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar,cutstar=cutstar))
    }
    else {
      return(list(mstar = mstar, sstar = sstar,cutstar=cutstar))
    }    
  } else {
    
    a   <- sum(p * mprime)
    b   <- sqrt(sum(p * ( (mprime - a)^2 + sprime ^2)))
    
    icchat_pw <- sum(p * ( (mprime - a)^2)) / b^2
    
    mstar <- (mprime - a) / b
    sstar <- sprime / b
    cutstar <- (cutests-a)/b 
    
    if (geticc == TRUE) {
      return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar,cutstar=cutstar))
    }
    else {
      return(list(mstar = mstar, sstar = sstar,cutstar=cutstar))
    } 
  }
  
}
