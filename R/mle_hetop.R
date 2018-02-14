###########################################################################################
## function to get MLE of (mug, sigmag, cutpoints) given first two fixed cutpoints.
## for unconstrained optimization parameters are organized as
## mug(G), log(sigmag)(G), log(cutpoint[3] - cutpoint[2), log(cutpoint[4] - cutpoint[3]) etc (K-3)
##
## updated 2/13/2018 to make a more complete accounting of data cases with
## respect to estimated parameters, to make sure MLE exists.
## Applies the following rules:
##
## 1) if group has >= 3 populated cells, estimate (mug, sigmag)
##
## 2) if a group has exactly 2 populated cells, estimate mug, and set log(sigmag)
##    = mean(log(sigmag)) of groups from case (1)
##
## 3) if a group has data in only 1 interior cell, proceed as with case (2)
##
## 4) if a group has data only in the top cell, set log(sigmag) = mean(log(sigmag))
##    of groups from case (1), and set mug = max of estimated means (cases 1+2+3)
##
## 5) if a group has data only in the bottom cell, set log(sigmag) = mean(log(sigmag))
##    of groups from case (1), and set mug = min of estimated means (cases 1+2+3)
###########################################################################################
mle_hetop <- function(ngk, fixedcuts, print.level=0, acov=FALSE, svals=NULL){
    
    stopifnot( is.matrix(ngk) && all(!is.na(ngk)) && all(ngk >= 0) && all(apply(ngk,1, sum) > 0) )
    G <- nrow(ngk)
    K <- ncol(ngk)
    if(K <= 2){
        stop("function currently required K >= 3 categories")
    }
    stopifnot( (length(fixedcuts) == 2) && (fixedcuts[1] < fixedcuts[2]) )

    ## define metadata "pstatus" giving, for each group, the estimation status of
    ## (mug, sigmag)
    pstatus <- as.data.frame(matrix("est", ncol=2, nrow=G), stringsAsFactors=FALSE)
    names(pstatus) <- c("mug","sigmag")
    
    pstatus$sigmag[which( apply(ngk, 1, function(x){ sum(x > 0) }) <= 2 )] <- "mean"
    pstatus$mug[which(apply(ngk, 1, function(x){ (sum(x > 0) == 1) && (x[1] > 0) }))] <- "min"
    pstatus$mug[which(apply(ngk, 1, function(x){ (sum(x > 0) == 1) && (x[K] > 0) }))] <- "max"

    n_m <- sum(pstatus$mug    == "est")
    n_s <- sum(pstatus$sigmag == "est")

    if(!( (n_m > 0) && (n_s > 0) )){
        stop("Insufficient data for estimation")
    }
    
    ## negative log likelihood, fixing parameters as needed
    negll_stz <- function(param){
        stopifnot(length(param) == (n_m + n_s + (K-3)))
        .mug <- .sigmag <- rep(-99.0, G)
        .mug[which(pstatus$mug == "est")] <- param[1:n_m]
        .mug[which(pstatus$mug == "min")] <- min(param[1:n_m])
        .mug[which(pstatus$mug == "max")] <- max(param[1:n_m])

        .lsigs <- param[(n_m + 1):(n_m + n_s)]
        .sigmag[which(pstatus$sigmag == "est")]  <- exp(.lsigs)
        .sigmag[which(pstatus$sigmag == "mean")] <- exp(mean(.lsigs))
        
        if(K==3){
            .cutpoints <- fixedcuts
        } else {
            .cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(n_m + n_s + 1):length(param)])))
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
        mu      <- mu[which(pstatus$mug == "est")]
        logsig  <- apply(ngk[which(pstatus$sigmag == "est"),], 1, function(x){ log(sqrt(weighted.mean(cats^2, w = x) - (weighted.mean(cats, w = x))^2)) })
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
    } else {
        acov <- NULL
    }
    
    if(K==3){
        cutpoints <- fixedcuts
    } else {
        cutpoints <- c(fixedcuts, fixedcuts[2] + cumsum(exp(param[(n_m + n_s + 1):length(param)])))
    }

    .mug <- .sigmag <- rep(-99.0, G)
    .mug[which(pstatus$mug == "est")] <- param[1:n_m]
    .mug[which(pstatus$mug == "min")] <- min(param[1:n_m])
    .mug[which(pstatus$mug == "max")] <- max(param[1:n_m])

    .lsigs <- param[(n_m + 1):(n_m + n_s)]
    .sigmag[which(pstatus$sigmag == "est")]  <- exp(.lsigs)
    .sigmag[which(pstatus$sigmag == "mean")] <- exp(mean(.lsigs))
    
    return(list(mug        = .mug,
                sigmag     = .sigmag,
                cutpoints  = cutpoints,
                acov       = acov,
                nlmdetails = res,
                pstatus    = pstatus))
}
