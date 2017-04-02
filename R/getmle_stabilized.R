###########################################################################################
## 9/2/2016: revised MLE function that uses sum to zero constraints (unweighted, and without
## switching group orders) and does two things to increase stability:
## 1) adds 0.5 to counts for groups that would/could cause MLE to not exist
## 2) uses better starting values
##
## for unconstrained optimization parameters are organized as mug (G-1),
## log(sigmag) (G-1), cutpoint[1], log(cutpoint[k]-cutpoint[k-1]) (K-2)
##
## for both sum to zero constraints, the first group is the holdout group
###########################################################################################
getmle_stabilized <- function(ngk, print.level=0, acov=FALSE, svals=NULL){
  G <- nrow(ngk)
  K <- ncol(ngk)

  ## determining locations of any cases with all data in extreme cells, which we use to adjust starting sigmas
  allextreme <- which(apply(ngk, 1, function(x){ (sum(x) == x[length(x)]) | (sum(x) == x[1]) }))
  
  ## stabilizing data
  bad <- which(apply(ngk, 1, function(x){ sum(x > 0) }) <= 2)
  if(length(bad) > 0){
      ngk[bad,] <- ngk[bad,] + 0.5
  }

  ## negative log likelihood function
  negll_stz <- function(param){
    stopifnot(length(param) == (K-1) + 2*(G-1))
    .mug       <- c(-sum(param[1:(G-1)]), param[1:(G-1)])
    .sigmag    <- exp(c(-sum(param[G:(2*G-2)]), param[G:(2*G-2)]))
    .cutpoints <- c(param[2*G-1], param[2*G-1]+ cumsum(exp(param[(2*G):length(param)])))
    
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

  ## set starting values based on means and SDs treating category labels as a scale
  ##
  ## NOTE: tried various things here with mixed success.  one issue is that if a group
  ## has all data in an extreme cell, the sample variance will be small but the MLE variance
  ## will be large, so we try to fix that by bumping up the starting variances for these cases  
  if(is.null(svals)){
      cats <- 1:K
      mu      <- apply(ngk, 1, function(x){ weighted.mean(cats, w = x) })
      logsig  <- apply(ngk, 1, function(x){ log(sqrt(weighted.mean(cats^2, w = x) - (weighted.mean(cats, w = x))^2)) })
      if(length(allextreme) > 0){
          logsig[allextreme] <- logsig[allextreme] + 2.0
      }
      ## logsig  <- rep(0, nrow(ngk))
      ## mu      <- 0.5 * (mu - mean(mu)) / sd(mu)
      mu          <- mu - mean(mu)
      mu[mu < -2] <- -2
      mu[mu > 2]  <-  2
      mu      <- (mu - mean(mu)) / sd(mu)
      logsig  <- logsig - mean(logsig)
      ## cuts based on marginal counts and assuming N(0,1)
      .cuts <- qnorm(cumsum(apply(ngk, 2, sum))/sum(ngk))[-K]
      svals   <- c(mu[-1], logsig[-1], .cuts[1], log(diff(.cuts)))
  }

  ## get MLE
  res <- nlm(f = negll_stz, hessian=acov, print.level=print.level, p = svals, iterlim=1500)
  
  if(res$code > 1){
    cat(paste("\n############# bad case #############\n"))
  }
  
  param <- res$estimate
  if(acov){
    acov <- solve(res$hessian)
    colnames(acov) <- rownames(acov) <- c(paste("mu",1:(G-1),sep=""), paste("logsigma",1:(G-1),sep=""), paste("cut",1:(K-1),sep=""))
  } else {
    acov <- NULL
  }
  return(list(mug        = c(-sum(param[1:(G-1)]), param[1:(G-1)]),
              sigmag     = exp(c(-sum(param[G:(2*G-2)]), param[G:(2*G-2)])),
              cutpoints  = c(param[2*G-1], param[2*G-1]+ cumsum(exp(param[(2*G):length(param)]))),
              acov       = acov,
              nlmdetails = res,
              badgroups  = bad,
              svals      = svals))
}

## #########################################################################
## ## checks on MLE functions
## #########################################################################
## 
## ##########################
## ## consistency
## ##########################
## set.seed(4022)
## G <- 4
## K <- 4
## Ng <- rep(1000000, G)
## theta <- list(mug       = c(1.0, 0.8, -0.8, -1.0),
##               sigmag    = exp(c(-1.3, 1.0, -0.5, 0.8)),
##               cutpoints = c(-1.2, 0.2, 1.2))
## 
## ngk <- gendata_hetop(G, K, Ng, theta$mug, theta$sigmag, theta$cutpoints)
## 
## thetahat <- getmle_stz(ngk)
## print(cbind(unlist(theta), unlist(thetahat[c("mug","sigmag","cutpoints")])))
## 
## thetahat <- getmle_stz_unweighted(ngk)
## print(cbind(unlist(theta), unlist(thetahat[c("mug","sigmag","cutpoints")])))
## 
## thetahat <- getmle_fixedcuts(ngk, fixedcuts = c(-1.2, 0.2))
## print(cbind(unlist(theta), unlist(thetahat[c("mug","sigmag","cutpoints")])))
## 
## thetahat <- getmle_stabilized(ngk)
## print(cbind(unlist(theta), unlist(thetahat[c("mug","sigmag","cutpoints")])))

################################################################################
## checking that it returns correct ordering even when it switches internally to
## put largest group first
################################################################################
## set.seed(4038)
## nsim <- 100
## res <- vector(nsim, mode="list")
## 
## for(s in 1:nsim){
##   G <- 10
##   K <- 4
##   Ng <- 20+round(100*runif(G), digits=0)
##   mug <- rnorm(G)
##   sigmag <- runif(G, 0.7, 1.4)
##   cutpoints <- c(-1, 0 , 1.5)
##   
##   ## generate data where MLE exists
##   ok <- FALSE
##   while(!ok){
##     ngk <- gendata_hetop(G, K, Ng, mug, sigmag, cutpoints)
##     if(min(apply(ngk, 1, function(x){ sum(x > 0)})) >= 2){
##       ok <- TRUE
##     }
##   }
##   
##   m0 <- getmle_stz(ngk, acov=TRUE, bigonefirst=FALSE, gradtol=1e-5)
##   m1 <- getmle_stz(ngk, acov=TRUE, bigonefirst=TRUE,  gradtol=1e-5)
##   res[[s]] <- data.frame(code0 = m0$nlmdetails$code,
##                          code1 = m0$nlmdetails$code,
##                          madiff = max(abs(unlist(m0[c("mug","sigmag","cutpoints")]) - unlist(m1[c("mug","sigmag","cutpoints")]))))
##   print(s)
## }
## 
## res <- do.call("rbind", res)
## print(sapply(res, summary))
## print(table(res$code0, res$code1))
## not clear this really improves anything

##
########################################
## checking parameterization invariance of MLE
########################################
## G <- 25
## K <- 4
## Ng <- c(1000, sample(20:200, size=G-1, replace=TRUE)) ## big diffs in weights
## theta <- list(mug       = rnorm(G, sd=0.7),
##               sigmag    = exp(rnorm(G, mean=0, sd=0.1)),
##               cutpoints = c(-1.5, 0.0, 1.5))
## 
## ngk <- gendata_hetop(G, K, Ng, theta$mug, theta$sigmag, theta$cutpoints)
## thetahat_u <- getmle_stz_unweighted(ngk, acov=FALSE)
## thetahat_w <- getmle_stz(ngk, acov=FALSE)
## p_g <- Ng /sum(Ng)
## print(a   <- sum(p_g * thetahat_u$mug))
## print(b   <- exp(sum(p_g * log(thetahat_u$sigmag))))
## 
## mug2 <- (thetahat_u$mug - a) / b
## print(summary(abs(thetahat_w$mug - mug2)))
## 
## sigmag2 <- thetahat_u$sigmag / b
## print(summary(abs(thetahat_w$sigmag - sigmag2)))
## 
## cuts <- (thetahat_u$cutpoints - a) / b
## print(summary(abs(thetahat_w$cutpoints - cuts)))
## 
## ## looks good up to numerical tolerance

