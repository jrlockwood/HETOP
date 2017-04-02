# try different approach to sum to zeros
###########################################################################################
## function to get MLE of (mug[-1], sigmag[-1], cutpoints)
## for unconstrained optimization parameters are organized as mug (G-1),
## log(sigmag) (G-1), cutpoint[1], log(cutpoint[k]-cutpoint[k-1]) (K-2)
##
## NOTE: here use weighted sum to zero constraints
###########################################################################################
getmle_stz_a <- function(ngk, print.level=0, acov=FALSE, svals=NULL){
  G <- nrow(ngk)
  K <- ncol(ngk)
  ng <- apply(ngk,1,sum) ##Kat added (2/28)
  pg <- ng/sum(ng) ##Kat added (2/28)
  
  negll_stz <- function(param){
    #stopifnot(length(param) == (K-1) + 2*(G-1))
    stopifnot(length(param) == (K-1) + 2*(G))
    ##    .mug       <- c(-sum(param[1:(G-1)]), param[1:(G-1)])
    ##    .sigmag    <- exp(c(-sum(param[G:(2*G-2)]), param[G:(2*G-2)]))
    ##Kat changed so that we constrain weighted sum = 0 instead of unweigthed sum
    # .mug       <- c(-sum(pg[2:G]*param[1:(G-1)])*(1/pg[1]), param[1:(G-1)])
    # .sigmag    <- exp(c(-sum(pg[2:G]*param[G:(2*G-2)])*(1/pg[1]), param[G:(2*G-2)]))
    # .cutpoints <- c(param[2*G-1], param[2*G-1]+ cumsum(exp(param[(2*G):length(param)])))
    
    .mug <- param[1:G]-(sum(pg * param[1:G]))
    .sigmag <- exp(c(param[(G+1):(2*G)] - sum(pg * param[(G+1):(2*G)])))
    .cutpoints <- c(param[2*G+1], param[2*G+1] + cumsum(exp(param[(2*G+2):length(param)])))
    
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
    #    svals <- c(rep(0,(G-1)),rep(0,(G-1)),seq(from=-1.5, to=1.5, length=K-1))
    svals <- c(rep(0,(G)),rep(0,(G)),seq(from=-1.5, to=1.5, length=K-1))
  }
  
  res <- nlm(f = negll_stz, hessian=acov, print.level=print.level, p = svals, iterlim=1500)
  
  if(res$code > 1){
    cat(paste("\n############# bad case #############\n"))
  }
  
  param <- res$estimate
  if(acov){
    acov <- solve(res$hessian)
    colnames(acov) <- rownames(acov) <- c(paste("mu",1:(G),sep=""), paste("logsigma",1:(G),sep=""), paste("cut",1:(K-1),sep=""))
  } else {
    acov <- NULL
  }
  ##Kat changed to constrain weighted sums (2/28)
  return(list( mug       = param[1:G]-(sum(pg * param[1:G])),
               sigmag    = exp(c(param[(G+1):(2*G)] - sum(pg * param[(G+1):(2*G)]))),
               cutpoints = c(param[2*G+1], param[2*G+1] + cumsum(exp(param[(2*G+2):length(param)]))),
               ##		  mug       = c(-sum(param[1:(G-1)]), param[1:(G-1)]),
               ##            sigmag    = exp(c(-sum(param[G:(2*G-2)]), param[G:(2*G-2)])),
               # cutpoints = c(param[2*G-1], param[2*G-1]+ cumsum(exp(param[(2*G):length(param)]))),
               acov      = acov,
               nlmdetails = res))
}
