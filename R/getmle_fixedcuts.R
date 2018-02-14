###########################################################################################
## function to get MLE of (mug, sigmag, cutpoints) given first two fixed cutpoints.
## for unconstrained optimization parameters are organized as
## mug(G), log(sigmag)(G), log(cutpoint[3] - cutpoint[2), log(cutpoint[4] - cutpoint[3]) etc (K-3)
##
## updated 05/02/2017 to set the log(sigamg) for groups with two or fewer populated
## cells to the mean(log(sigmag)) across other groups to force model identification
###########################################################################################
getmle_fixedcuts <- function(ngk, fixedcuts, print.level=0, acov=FALSE, svals=NULL){
  G <- nrow(ngk)
  K <- ncol(ngk)
  stopifnot( (length(fixedcuts) == 2) & (fixedcuts[1] < fixedcuts[2]) )

  ## locate "bad" groups with degenerate cells
  bad    <- which(apply(ngk, 1, function(x){ sum(x > 0) }) <= 2)
  notbad <- setdiff(1:G, bad)
  nbad   <- length(bad)
  if(nbad == G){
      stop("There are no groups with sufficient data to estimate sigma_g")
  }

  ## negative log likelihood, fixing sigmag = exp(mean(log(sigmag))) for bad groups
  negll_stz <- function(param){
    stopifnot(length(param) == G + (G-nbad) + (K-3) )
    .mug            <- param[1:G]
    .lsigs          <- param[(G+1):(2*G - nbad)]
    .sigmag         <- rep(-1.0,G)
    .sigmag[notbad] <- exp(.lsigs)
    .sigmag[bad]    <- exp(mean(.lsigs))

    if(K==3){
        .cutpoints <- fixedcuts
    } else {
        .cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(2*G - nbad + 1):length(param)])))
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
      logsig  <- apply(ngk[notbad,], 1, function(x){ log(sqrt(weighted.mean(cats^2, w = x) - (weighted.mean(cats, w = x))^2)) })
      mu      <- (mu - mean(mu)) / sd(mu)
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
      ## if(K==3){
      ##     colnames(acov) <- rownames(acov) <- c(paste("mu",1:G,sep=""), paste("logsigma",1:G,sep=""))
      ## } else {
      ##    colnames(acov) <- rownames(acov) <- c(paste("mu",1:G,sep=""), paste("logsigma",1:G,sep=""), paste("cutdev",1:(K-3),sep=""))
      ## }
  } else {
    acov <- NULL
  }

  if(K==3){
      cutpoints <- fixedcuts
  } else {
      cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(2*G - nbad + 1):length(param)])))
  }

  .lsigs          <- param[(G+1):(2*G - nbad)]
  .sigmag         <- rep(-1.0,G)
  .sigmag[notbad] <- exp(.lsigs)
  .sigmag[bad]    <- exp(mean(.lsigs))

  return(list(mug       = param[1:G],
              sigmag    = .sigmag,
              cutpoints = cutpoints,
              acov      = acov,
              nlmdetails = res))
}
