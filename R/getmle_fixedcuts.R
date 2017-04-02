###########################################################################################
## function to get MLE of (mug, sigmag, cutpoints) given first two fixed cutpoints.
## for unconstrained optimization parameters are organized as
## mug(G), log(sigmag)(G), log(cutpoint[3] - cutpoint[2), log(cutpoint[4] - cutpoint[3]) etc (K-3)
###########################################################################################
getmle_fixedcuts <- function(ngk, fixedcuts, print.level=0, acov=FALSE, svals=NULL){
  G <- nrow(ngk)
  K <- ncol(ngk)
  stopifnot( (length(fixedcuts) == 2) & (fixedcuts[1] < fixedcuts[2]) )
  
  negll_stz <- function(param){
    stopifnot(length(param) == 2*G + (K-3))
    .mug       <- param[1:G]
    .sigmag    <- exp(param[(G+1):(2*G)])
    if(K==3){
        .cutpoints <- fixedcuts
    } else {
        .cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(2*G+1):length(param)])))
    }
    
    ## calculate cell probabilities, truncating at values very close to 0 or 1
    pgk <- matrix(0, ncol=K, nrow=G)
    for(g in 1:G){
      tmp     <- c(pnorm(.cutpoints, mean=.mug[g], sd = .sigmag[g]), 1)
      pgk[g,] <- c(tmp[1], diff(tmp))
    }
    stopifnot(max(abs(apply(pgk, 1, sum) - 1)) < 1e-5)
    pgk[which(pgk <= 1e-300)]     <- 1e-300
    pgk[which(pgk >= 1 - 1e-300)] <- 1 - 1e-300
    
    ## return negative of multinomial log likelihood based on counts "ngk"
    -sum(c(ngk)*log(c(pgk)))
  }

  if(is.null(svals)){
      cats <- 1:K
      mu      <- apply(ngk, 1, function(x){ weighted.mean(cats, w = x) })
      logsig  <- apply(ngk, 1, function(x){ log(sqrt(weighted.mean(cats^2, w = x) - (weighted.mean(cats, w = x))^2)) })
      mu      <- (mu - mean(mu)) / sd(mu)
      logsig[which(logsig < -4)] <- -4 ## groups with degenerate counts
      logsig  <- logsig - mean(logsig)
      svals   <- c(mu, logsig, rep(0, K-3))
  }
  
  res <- nlm(f = negll_stz, hessian=acov, print.level=print.level, p = svals, iterlim=1500)
  
  if(res$code > 1){
    cat(paste("\n############# bad case #############\n"))
  }
  
  param <- res$estimate
  if(acov){
      acov <- solve(res$hessian)
      if(K==3){
          colnames(acov) <- rownames(acov) <- c(paste("mu",1:G,sep=""), paste("logsigma",1:G,sep=""))
      } else {
          colnames(acov) <- rownames(acov) <- c(paste("mu",1:G,sep=""), paste("logsigma",1:G,sep=""), paste("cutdev",1:(K-3),sep=""))
      }
  } else {
    acov <- NULL
  }

  if(K==3){
      cutpoints <- fixedcuts
  } else {
      cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(2*G+1):length(param)])))
  }
  
  return(list(mug       = param[1:G],
              sigmag    = exp(param[(G+1):(2*G)]),
              cutpoints = cutpoints,
              acov      = acov,
              nlmdetails = res))
}
