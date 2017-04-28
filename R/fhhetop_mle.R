## (4/27/2017): MLE estimation of FH-HETOP model

## NOTE: THIS VERSION REQUIRES at least one covariate for means and logSDs to
## avoid too many code forks at this development stage

## NOTE: using large M is not computationally feasible with this code

## NOTE: As Efron notes, MLEs for alpha are ill-behaved; he uses penalty, here
## we range restrict using "alphamin" and "alphamax"

fhhetop_mle <- function(ngk, Xm=NULL, Xs=NULL, fixcuts = c(-1.0, 0.0), p, m, gridL, gridU, print.level=0, acov=FALSE, svals=NULL, alphamin=-50, alphamax = 50){

    expit <- function(x){
        v<-exp(x)
        v/(1+v)
    }

    stopifnot(is.matrix(ngk) && all(!is.na(ngk)) && all(c(ngk) >= 0))
    G <- nrow(ngk)
    K <- ncol(ngk)
    Kminus3 <- K-3
  
    if(K < 3){
        stop("current implementation requires at least 3 categories")
    }
        
    ng	<- apply(ngk, 1, sum)
    pg	<- ng/sum(ng)

    dimXm <- dimXs <- 0
    meanX <- !is.null(Xm)  
    if(meanX){
        stopifnot(is.matrix(Xm) && all(!is.na(Xm)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
        dimXm <- ncol(Xm)
    }
    sdX   <- !is.null(Xs)
    if(sdX){
        stopifnot(is.matrix(Xs) && all(!is.na(Xs)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
        dimXs <- ncol(Xs)
    }
    
    ## if(as.integer(meanX + sdX) == 1){
    ##    stop("current implementation requires mean structure be specified for neither or both of (mu, sigma)")
    ## }
    if(!meanX | !sdX){
        stop("this devel version requires at least one covariate for each of means and log(sds)")
    }
    
    cuts12 <- sort(fixcuts)
    stopifnot( (length(cuts12)==2) && all(!is.na(cuts12)) && is.numeric(cuts12) )
    
    ############################
    ## set up spline bases and bivariate grid
    ############################
    gridm  <- seq(from = gridL[1], to = gridU[1], length=m[1])
    Qm     <- bs(gridm, df = p[1], degree=3, intercept=FALSE)

    grids  <- seq(from = gridL[2], to = gridU[2], length=m[2])
    Qs     <- bs(grids, df = p[2], degree=3, intercept=FALSE)

    gridms <- expand.grid(gridm, grids)

    #############################
    ## negative log likelihood for cuts, regression coeffs, Efron prior parameters
    #############################
    negll <- function(param){
        
        if(K==3){
            cuts <- cuts12
        } else {
            cuts <- c(cuts12, cumsum(exp(param[1:Kminus3])))
        }
        beta_m  <- param[(Kminus3 + 1):(Kminus3 + dimXm)]
        beta_s  <- param[(Kminus3 + dimXm + 1):(Kminus3 + dimXm + dimXs)]
        alpha_m <- alphamin + expit(param[(Kminus3 + dimXm + dimXs + 1):(Kminus3 + dimXm + dimXs + p[1])])*(alphamax - alphamin)
        alpha_s <- alphamin + expit(param[(Kminus3 + dimXm + dimXs + p[1] + 1):(Kminus3 + dimXm + dimXs + p[1] + p[2])])*(alphamax - alphamin)
        .gamma  <- param[(Kminus3 + dimXm + dimXs + p[1] + p[2] + 1)]

        ## get probabilities on bivariate grid
        pFm     <- exp(Qm %*% alpha_m)
        pFm     <- pFm / sum(pFm)        
        pFs     <- exp(Qs %*% alpha_s)
        pFs     <- pFs / sum(pFs)
        pF      <- apply(expand.grid(pFm, pFs),1,prod)
        
        ## sum log integrated likelihood over rows of ngk
        ## NOTE: uses trick for dealing with very small numbers
        ll <- sapply(1:G, function(g){
            sigmag <- exp(gridms[,2] + sum(Xs[g,] * beta_s)) ## note: this is redundant, could be improved
            mug    <- gridms[,1] + .gamma*gridms[,2] + sum(Xm[g,] * beta_m)
            llcond <- sapply(1:nrow(gridms), function(i){
                tmp  <- c(pnorm(cuts, mean=mug[i], sd = sigmag[i]), 1)
                pgk  <- c(tmp[1], diff(tmp))
                pgk[which(pgk <= 1e-300)]     <- 1e-300
                pgk[which(pgk >= 1 - 1e-300)] <- 1 - 1e-300
                sum(ngk[g,] * log(pgk))
            })
            .max <- max(llcond)
            return( .max + log(sum(pF * exp(llcond - .max))) )
        })
        return(-sum(ll))
    }

    ############################
    ## set starting values and do optimization
    ############################
    if(is.null(svals)){
        cats <- 1:K
        mu      <- apply(ngk, 1, function(x){ weighted.mean(cats, w = x) })
        logsig  <- apply(ngk, 1, function(x){ log(sqrt(weighted.mean(cats^2, w = x) - (weighted.mean(cats, w = x))^2)) })
        logsig[!is.finite(logsig)] <- mean(logsig[is.finite(logsig)])
        mu      <- as.vector(scale(mu))

        svals    <- c(rep(log(1),Kminus3),
                      as.vector(coef(lm(mu ~ Xm - 1))),
                      as.vector(coef(lm(logsig ~ Xs -1))),
                      rep(0.0, p[1]),
                      rep(0.0, p[2]),
                      0.0)
    }

    sol <- nlm(f = negll, hessian=acov, print.level=print.level, p = svals, iterlim=1500)

    #####################################
    ## return pieces
    #####################################
    param <- sol$estimate
    return(list(nlminfo = sol,
                beta_m  = param[(Kminus3 + 1):(Kminus3 + dimXm)],
                beta_s  = param[(Kminus3 + dimXm + 1):(Kminus3 + dimXm + dimXs)],
                alpha_m = alphamin + expit(param[(Kminus3 + dimXm + dimXs + 1):(Kminus3 + dimXm + dimXs + p[1])])*(alphamax - alphamin),
                alpha_s = alphamin + expit(param[(Kminus3 + dimXm + dimXs + p[1] + 1):(Kminus3 + dimXm + dimXs + p[1] + p[2])])*(alphamax - alphamin),
                .gamma  = param[(Kminus3 + dimXm + dimXs + p[1] + p[2] + 1)]))
}



## LOAD SIM DATA TO GET A TEST CASE
## load("~/Projects/2015-HETOP-model/fay-herriot-code/performance-sims/from-cluster/efron/N100-skewed-messy-p1010.Robj")
## ex      <- simres[[1]]$Dinfo
## 
## ngk     <- ex$ngk
## Xm      <- Xs <- matrix(ex$xg - mean(ex$xg), ncol=1)
## fixcuts <- c(-1,0)
## p       <- c(5,5)
## m       <- c(20,20)
## gridL   <- c(-5.0, log(0.10))
## gridU   <- c( 5.0, log(5.00))
## 
## tmp     <- fhhetop_mle(ngk, Xm, Xs, fixcuts, p, m, gridL, gridU, print.level=2, acov=FALSE, svals=NULL)
