############################################################
## (7/12/2016): FH-HETOP model with Efron prior
## (updated in v19 to use different alpha parameterization that seems to mix better, drop to 2 chains)
############################################################
hetop_fh_jags_efron_v1 <- function(ngk, Xm=NULL, Xs=NULL, bugs.seeds, nburn, niter, useglm=FALSE, fixcuts = c(-1.0, 0.0), workdir=NULL, cleandir=FALSE, p, m, gridL, gridU){
    ## inherits from hetop_s1_fh_jags_v2f() but uses efron prior and linear regression on residuals.
    ## additional arguments needed here are.
    ##
    ## p:     2-dimensional vector of degrees of freedom for efron prior
    ## m:     2-dimensional vector of number of grid points for efron prior
    ## gridL: 2-dimensional vector of lower boundaries for grid
    ## gridU: 2-dimensional vector of upper boundaries for grid
  
    stopifnot(is.matrix(ngk) && all(!is.na(ngk)) && all(c(ngk) >= 0))
  
    if(ncol(ngk) < 3){
        stop("current implementation requires at least 3 categories")
    }
        
    ng	<- apply(ngk, 1, sum)
    pg	<- ng/sum(ng)
    
    meanX <- !is.null(Xm)  
    if(meanX){
        stopifnot(is.matrix(Xm) && all(!is.na(Xm)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
    }
    sdX   <- !is.null(Xs)
    if(sdX){
        stopifnot(is.matrix(Xs) && all(!is.na(Xs)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
    }
    
    if(as.integer(meanX + sdX) == 1){
        stop("current implementation requires mean structure be specified for neither or both of (mu, sigma)")
    }
    
    cuts12 <- sort(fixcuts)
    stopifnot( (length(cuts12)==2) && all(!is.na(cuts12)) && is.numeric(cuts12) )
    
    if (is.null(workdir)) {
        stop("need to specify a directory for storing model files and results")
    } else {
        setwd(workdir)
        if (cleandir) {
            for (fn in c("jagschain0", "jagschain1", "jagschain2", "jagschain3", "jagsdata",
                         "jagsindex","jagsindex0", "jagsindex1", "jagsindex2", "jagsindex3",
                         "jagsinits1", "jagsinits2", "jagsinits3",
                         "jagsmodel", "jagsscript")) {
                if(file.exists(paste0(fn, ".txt"))) file.remove(paste0(fn, ".txt"))
            }
        }
    }

    ############################
    ## set up spline bases (need to strip attributes from Qm and Qs)
    ############################
    gridm  <- seq(from = gridL[1], to = gridU[1], length=m[1])
    Qm     <- bs(gridm, df = p[1], degree=3, intercept=FALSE)
    Qm     <- matrix(c(Qm), ncol=p[1], byrow=F)

    grids  <- seq(from = gridL[2], to = gridU[2], length=m[2])
    Qs     <- bs(grids, df = p[2], degree=3, intercept=FALSE)
    Qs     <- matrix(c(Qs), ncol=p[2], byrow=F)
    
    #####################
    ## create data file
    #####################
    dat <- list(G = nrow(ngk), K = ncol(ngk), ngk = ngk, N = apply(ngk, 1, sum), m = m, p = p, Qm = Qm, Qs = Qs, gridm = gridm, grids = grids)
    if(meanX){
        dat$Xm    <- Xm
        dat$dimXm <- ncol(Xm)
    }
    if(sdX){
        dat$Xs    <- Xs
        dat$dimXs <- ncol(Xs)
    }
    if(dat$K <= 3){
        dat$cuts <- cuts12
    }
    if(dat$K >= 4){
        dat$cuts <- c(cuts12, rep(NA, dat$K-3))
    }
    dump(names(dat), envir=as.environment(dat), file = "jagsdata.txt")
    
    #####################
    ## create model file
    ## NOTE: some forks necessary
    #####################
    cat("model
  {
    for(g in 1:G){
      ngk[g,1:K] ~ dmulti(pgk[g,1:K], N[g])

      pgk[g,1] <- phi( (cuts[1] - mu[g]) / sigma[g])
      for(k in 2:(K-1)){
        pgk[g,k] <- phi( (cuts[k] - mu[g]) / sigma[g]) - sum(pgk[g,1:(k-1)])
      }
      pgk[g,K] <- 1.0 - phi( (cuts[K-1] - mu[g]) / sigma[g])

      mu[g]    <-     epsilon[g,1]
      sigma[g] <- exp(epsilon[g,2])

      locm[g] ~ dcat(pFm[1:m[1]])
      locs[g] ~ dcat(pFs[1:m[2]])\n", file="jagsmodel.txt")
    
    if( (!meanX && !sdX) ){
        cat("
      epsilon[g,1]  <- (gamma * grids[locs[g]]) + gridm[locm[g]]
      epsilon[g,2]  <- grids[locs[g]]\n", file="jagsmodel.txt", append=TRUE)
    }
    if( (meanX && sdX) ){
        cat("
      epsilon[g,1] <- inprod(Xm[g,1:dimXm], beta_m[1:dimXm]) + (gamma * grids[locs[g]]) + gridm[locm[g]]
      epsilon[g,2] <- inprod(Xs[g,1:dimXs], beta_s[1:dimXs]) + grids[locs[g]]\n", file="jagsmodel.txt", append=TRUE)
    }
    cat("    }\n", file="jagsmodel.txt", append=TRUE)

    cat("
    alpha0m     ~ dnorm(0.0, 0.01)
    for(i in 1:p[1]){
       alphamdev[i] ~ dnorm(0.0, 0.05)
       alpham[i] <- alpha0m + alphamdev[i]
    }
    for(i in 1:m[1]){
      pFm[i] <- exp(Qm[i,1:p[1]] %*% alpham[1:p[1]])
    }

    alpha0s     ~ dnorm(0.0, 0.01)
    for(i in 1:p[2]){
       alphasdev[i] ~ dnorm(0.0, 0.05)
       alphas[i] <- alpha0s + alphasdev[i]
    }
    for(i in 1:m[2]){
      pFs[i] <- exp(Qs[i,1:p[2]] %*% alphas[1:p[2]])
    }

    gamma ~ dnorm(0.0, 0.1)

    ", file="jagsmodel.txt", append=TRUE) 

    if(dat$K >= 4){
        cat("
    for(k in 1:(K-3)){
      cuts0[k] ~ dnorm(0.0, 0.01) I(cuts[2], )
    }
    tmp[1:(K-3)] <- sort(cuts0)
    for(k in 1:(K-3)){
      cuts[k+2] <- tmp[k]
    }
    ", file="jagsmodel.txt", append=TRUE)
    }

    if(meanX){
        cat("
    for(i in 1:dimXm){
      beta_m[i] ~ dnorm(0.0, 0.1)
    }
    ", file="jagsmodel.txt", append=TRUE)
    }
    if(sdX){
        cat("
    for(i in 1:dimXs){
      beta_s[i] ~ dnorm(0.0, 0.1)
    }
    ", file="jagsmodel.txt", append=TRUE)
    }
    cat("}\n", file="jagsmodel.txt", append=TRUE)
    
    #####################
    ## generate inits
    ##
    ## NOTE: ran into problems with JAGS complaining about data being inconsistent
    ## with initial params - think this is due to starting SDs being small which
    ## provide nearly zero probabilities for some cells that actually have data.
    ## so fix this by starting SDs at 80th percentile.  also disperse starting values
    ## more by starting means at random locations.
    ##
    ## NOTE: found this problem creeping up again even starting at the 80th percentile
    ## so we switch to putting them at the 90th percentile
    #####################
    gen.inits <- function(i, seed){
        .tmp <- list(locm      = as.integer(sample(1:m[1], size=dat$G, replace=TRUE)),
                     locs      = as.integer(rep(floor(0.90 * m[2]), dat$G)),
                     alpha0m   = rnorm(1),
                     alphamdev = rnorm(p[1]),
                     alpha0s   = rnorm(1),
                     alphasdev = rnorm(p[2]),
                     gamma     = 0.0,
                     .RNG.name = "base::Mersenne-Twister",
                     .RNG.seed = seed)
        if(meanX){
            .tmp$beta_m <- rnorm(dat$dimXm, sd=0.1)
        }
        if(sdX){
            .tmp$beta_s <- rnorm(dat$dimXs, sd=0.1)
        }
        if(dat$K >= 4){
            .tmp$cuts0 <- sort(seq(from = cuts12[2] + 1, to = cuts12[2] + 3, length=dat$K-3))
        }
        
        dump(names(.tmp), envir=as.environment(.tmp), file=paste("jagsinits",i,".txt",sep=""))
    }
  
    gen.inits(1, bugs.seeds[1])
    gen.inits(2, bugs.seeds[2])
    ## gen.inits(3, bugs.seeds[3])

    ####################
    ## write script
    ####################
    cat("load dic\n", file="jagsscript.txt")
    if(useglm){
        cat("load glm\n", file="jagsscript.txt", append=TRUE)
    }
    cat("model in 'jagsmodel.txt'
data in 'jagsdata.txt'
compile, nchains(2)
parameters in 'jagsinits1.txt', chain(1)
parameters in 'jagsinits2.txt', chain(2)
initialize\n", file="jagsscript.txt", append=TRUE)
    cat(paste("update ",nburn,", by (20)\n",sep=""), file="jagsscript.txt", append=TRUE)
    cat("monitor set mu
monitor set sigma
monitor set cuts\n", file="jagsscript.txt", append=TRUE)
    if(meanX){
        cat("monitor set beta_m\n", file="jagsscript.txt", append=TRUE)
    }
    if(sdX){
        cat("monitor set beta_s\n", file="jagsscript.txt", append=TRUE)
    }
    cat("monitor set alpha0m
monitor set alpham
monitor set alpha0s
monitor set alphas
monitor set gamma
monitor set pD
monitor set deviance\n", file="jagsscript.txt", append=TRUE)
    cat(paste("update ",niter,", by (20)\n",sep=""), file="jagsscript.txt", append=TRUE)
    cat("coda *, stem(jags)\n", file = "jagsscript.txt", append=TRUE)
  
    ############################
    ## run JAGS and get results
    ############################
    system("jags jagsscript.txt")
    
    ## note: pD information is stored in its own CODA file because it pools across
    ## chains
    pD <- read.coda("jagschain0.txt","jagsindex0.txt", quiet=TRUE)
    r1 <- read.coda("jagschain1.txt", "jagsindex.txt", quiet=TRUE)
    r2 <- read.coda("jagschain2.txt", "jagsindex.txt", quiet=TRUE)
    ## r3 <- read.coda("jagschain3.txt", "jagsindex.txt", quiet=TRUE)
    ## rall <- rbind(r1, r2, r3)
    rall <- rbind(r1, r2)
    
    ## note deviance is at end, grab it and take it off and do DIC calc
    ## DIC = mean(deviance) + mean(pD)
    ## NOTE: from thread here:
    ## http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
    deviance <- rall[,"deviance"]
    rall <- rall[,setdiff(colnames(rall),"deviance")]
    dic <- mean(deviance) + mean(pD)
    
    gr <- vector(ncol(rall), mode="list")
    for(param in 1:ncol(rall)){
        ## gr[[param]] <- gelman.diag(mcmc.list(r1[,param],r2[,param], r3[,param]))$psrf
        gr[[param]] <- gelman.diag(mcmc.list(r1[,param],r2[,param]))$psrf
    }
    gr <- as.data.frame(do.call("rbind", gr))
    names(gr) <- c("gr.est","gr.upper")
    r <- data.frame(mean = apply(rall, 2, mean), sd = apply(rall, 2, sd), gr)
    cat("\n\nPARAMETERS WITH HIGH GR\n")
    print(subset(r, gr.est >= 1.05), digits=4)
    
    return(list(samps = rall, summaries = r, deviance = deviance, pD = pD, dic = dic, Qm = Qm, Qs = Qs))
}
