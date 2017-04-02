##########################################################################################
## function to compute WAIC from JAGS samples
## formula taken from Gelman et al paper
##########################################################################################
waic_hetop <- function(ngk, samps){
  stopifnot(is.matrix(ngk) && all(!is.na(ngk)))  
  G <- nrow(ngk)
  K <- ncol(ngk)

  ## get blocks of parameter samples
  mug <- samps[,grep("mu",colnames(samps))]
  stopifnot(ncol(mug) == G)
  sigmag <- samps[,grep("sigma",colnames(samps))]
  stopifnot(ncol(sigmag) == G)
  cuts <- samps[,grep("cuts",colnames(samps))]
  stopifnot(ncol(cuts) == (K-1))

  ## compute (nsim x G) matrix of log likelihoods
  ll <- t(sapply(1:nrow(samps), function(s){
      pgk <- matrix(0, ncol=K, nrow=G)
      for(g in 1:G){
          tmp     <- c(pnorm(cuts[s,], mean=mug[s,g], sd = sigmag[s,g]), 1)
          pgk[g,] <- c(tmp[1], diff(tmp))
      }
      pgk[which(pgk <= 1e-300)]     <- 1e-300
      apply(ngk * log(pgk), 1, sum)
  }))

  lpd_hat <- sum(log(apply(ll, 2, function(logparts){
      m <- max(logparts)
      exp(m) * mean(exp(logparts - m))
  })))

  phat_waic <- sum(apply(ll, 2, var))

  return(list(lpd_hat = lpd_hat, phat_waic = phat_waic, waic = -2 * (lpd_hat - phat_waic)))
}
