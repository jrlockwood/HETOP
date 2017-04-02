##KAT ADDED
##3/20/2016: KAT MODIFIED TO PUT REG COEFFS ON STAR METRIC
##8/14/2016: JRL added return_ab argument
est_deref_jags <- function(mprime, sprime, cutests, ngs, geticc, getregs=FALSE, beta_ms=NULL, beta_ss=NULL, return_ab=FALSE) {
    ##   function to "standardize" the JAGS estimates
    ##   mus=vector of mean estimates
    ##   sigs=vector of sigma estimates
    ##   ngs=vector of group sample sizes
    ##   getregs=logical to indicate if reg coeffs should be standardized and returned
    ## 	 beta_ms = vector of reg coeffs for mns
    ##   beta_ss = vector of reg coefs for sds 
    
    N <- sum(ngs)
    G <- length(ngs)
    p <- ngs/N
    
    a   <- sum(p * mprime)
    b   <- sqrt(sum(p * ( (mprime - a)^2 + sprime ^2)))
    
    icchat_pw <- sum(p * ( (mprime - a)^2)) / b^2
    
    mstar   <- (mprime - a) / b
    sstar   <- sprime / b
    cutstar <- (cutests-a)/b  
    
    if(getregs == TRUE) {
        b_mstar	<- beta_ms/b
        b_sstar	<- beta_ss
        if(geticc == TRUE) {
            if(return_ab){
                return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar,cutstar = cutstar, b_mstar = b_mstar, b_sstar = b_sstar, a = a, b = b))
            } else {
                return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar,cutstar = cutstar, b_mstar = b_mstar, b_sstar = b_sstar))
            }
        } else {
            if(return_ab){
                return(list(mstar = mstar, sstar = sstar,cutstar = cutstar, b_mstar = b_mstar, b_sstar = b_sstar, a = a, b = b))
            } else {
                return(list(mstar = mstar, sstar = sstar,cutstar = cutstar, b_mstar = b_mstar, b_sstar = b_sstar))
            }
        }
    }
    if(getregs==FALSE) {
        if (geticc == TRUE) {
            if(return_ab){
                return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar, cutstar=cutstar, a = a, b = b))
            } else {
                return(list(icchat_pw = icchat_pw, mstar = mstar, sstar = sstar, cutstar=cutstar))
            }
        }
        else {
            if(return_ab){
                return(list(mstar = mstar, sstar = sstar, cutstar = cutstar, a = a, b = b))
            } else {
                return(list(mstar = mstar, sstar = sstar, cutstar = cutstar))
            }
        }
    }
}
