#####################################################################
## (4/6/16): JAGS function for better mixing of latent regression
## (4/11/16): add options to control the working directory it stores files in and cleanup directory
## (7/12/16): add option to fix two cutpoints and estimate means (default FALSE)
#####################################################################
hetop_s1_fh_jags_v2f <- function(ngk, Xm=NULL, Xs=NULL, bugs.seeds, nburn, niter, useglm=FALSE, fixcuts=NULL, svals=NULL, workdir=NULL, cleandir=FALSE){
    ################################################################
    ## Function to estimate single-subject Fay-Herriot HETOP model in JAGS.
    ## Replaces previous functions by covering both the no-covariate and covariate
    ## case, also allowing mean and log(SD) to differ on whether covariates are
    ## included.  So this really covers four cases.
    ##
    ## Also fixes the hack of combining sum to zero with random effects because I
    ## figured out that provided non-exchangeable priors so there would be label
    ## dependence. This version does not force sum to zero but encourages it by
    ## using a LR model that has a true zero intercept because the X's are centered,
    ## and using a residual distribution with mean zero
    ##
    ## ngk:        G x K matrix of counts
    ## Xm:         G x dimXm matrix of covariates for means (if NULL, no covs)
    ## Xs:         G x dimXs matrix of covariates for log(sigmas) (if NULL, no covs)
    ## bugs.seeds: vector of at least three integer seeds for set.seed and JAGS
    ## nburn:      number of burn-in iters per chain (3)
    ## niter:      number of saved iters per chain (3)
    ## useglm:     if TRUE, load glm module
    ## fixcuts:    if !NULL, uses alt identification of fixing 2 cutpoints and allowing nonzero means
    ##             in this case, fixcuts must a 2-dimensional vector of the leftmost two cutpoints
    ## svals:      starting values for (mug, sigmag, cuts) passed as list (if NULL, generate vals)
    ##             only if !fixcuts   
    ## workdir:    directory where model files and output will get stored
    ## cleandir:   if TRUE, this will try to erase the various files that would get created in the workdir 
    ##             directory during function call
    ###################################################################
    stopifnot(is.matrix(ngk) && all(!is.na(ngk)) && all(c(ngk) >= 0))
  
    if(ncol(ngk) < 3){
        stop("current implementation requires at least 3 categories")
    }
    
    if(!is.null(fixcuts)){
        cuts12 <- sort(fixcuts)
        stopifnot(length(cuts12)==2)
        fixcuts <- TRUE
    } else {
        fixcuts <- FALSE
    }
    
    if(!is.null(svals)){
        stopifnot(!fixcuts & is.list(svals) && all(names(svals)==c("mug","sigmag","cutpoints")))
    }
    
    ng	<- apply(ngk, 1, sum)
    pg	<- ng/sum(ng)
    
    ## Kat changed so checks that weighted means = 0 instead of unweighted means (2/28)
    ## JRL reversed this change 7/12/16
    meanX <- !is.null(Xm)  
    if(meanX){
        stopifnot(is.matrix(Xm) && all(!is.na(Xm)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
        ## stopifnot(is.matrix(Xm) && all(!is.na(Xm)) && max(abs(apply(Xm, 2, function(x) sum(pg*x)))) < 1e-6)
    }
    sdX   <- !is.null(Xs)
    if(sdX){
        stopifnot(is.matrix(Xs) && all(!is.na(Xs)) && max(abs(apply(Xs, 2, mean))) < 1e-6)
        ## stopifnot(is.matrix(Xs) && all(!is.na(Xs)) && max(abs(apply(Xm, 2, function(x) sum(pg*x)))) < 1e-6)
    }
    
    if(as.integer(meanX + sdX) == 1){
        stop("current implementation requires mean structure be specified for neither or both of (mu, sigma)")
    }
    
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
  
    #####################
    ## create data file
    #####################
    dat <- list(G = nrow(ngk), K = ncol(ngk), ngk = ngk, N = apply(ngk, 1, sum))
    if(meanX){
        dat$Xm    <- Xm
        dat$dimXm <- ncol(Xm)
    }
    if(sdX){
        dat$Xs    <- Xs
        dat$dimXs <- ncol(Xs)
    }
    if( (!meanX && !sdX) && !fixcuts ){
        dat$zeros <- c(0,0)
    }
    if(fixcuts && dat$K <= 3){
        dat$cuts <- cuts12
    }
    if(fixcuts && dat$K >= 4){
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
      sigma[g] <- exp(epsilon[g,2])\n", file="jagsmodel.txt")

    ## fork on four cases
    if( (!meanX && !sdX) & !fixcuts){
        cat("
      epsilon[g,1:2] ~ dmnorm(zeros[1:2], Sigmainv[1:2,1:2])\n", file="jagsmodel.txt", append=TRUE)
    }
    if( (!meanX && !sdX) & fixcuts){
        cat("
      epsilon[g,1:2] ~ dmnorm(mu0[1:2], Sigmainv[1:2,1:2])\n", file="jagsmodel.txt", append=TRUE)
    }
    if( (meanX && sdX) & !fixcuts){
        cat("
      meanstruct[g,1] <- inprod(Xm[g,1:dimXm], beta_m[1:dimXm])
      meanstruct[g,2] <- inprod(Xs[g,1:dimXs], beta_s[1:dimXs])
      epsilon[g,1:2] ~ dmnorm(meanstruct[g,1:2], Sigmainv[1:2,1:2])\n", file="jagsmodel.txt", append=TRUE)
    }
    if( (meanX && sdX) & fixcuts){
        cat("
      meanstruct[g,1] <- mu0[1] + inprod(Xm[g,1:dimXm], beta_m[1:dimXm])
      meanstruct[g,2] <- mu0[2] + inprod(Xs[g,1:dimXs], beta_s[1:dimXs])
      epsilon[g,1:2] ~ dmnorm(meanstruct[g,1:2], Sigmainv[1:2,1:2])\n", file="jagsmodel.txt", append=TRUE)
    }
    cat("    }\n", file="jagsmodel.txt", append=TRUE)
    
    ## deal with cuts depending on fixcuts; fixcuts==TRUE case is messy because the code must depend on K
    if( !fixcuts ){
        cat("
    for(k in 1:(K-1)){
      cuts0[k] ~ dnorm(0.0, 0.01)
    }
    cuts[1:(K-1)] <- sort(cuts0)
   ", file="jagsmodel.txt", append=TRUE)
    }

    if( fixcuts ){
        cat("
    mu0[1]  ~ dnorm(0.0, 0.01)
    mu0[2]  ~ dnorm(0.0, 0.01)
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
    cat("
    eps_sd[1] ~ dunif(0,2.0)
    eps_sd[2] ~ dunif(0,1.0)
    eps_corr  ~ dunif(-1,1)
    Sigma[1,1] <- eps_sd[1] * eps_sd[1]
    Sigma[1,2] <- eps_sd[1] * eps_sd[2] * eps_corr
    Sigma[2,1] <- eps_sd[1] * eps_sd[2] * eps_corr    
    Sigma[2,2] <- eps_sd[2] * eps_sd[2]
    Sigmainv[1:2,1:2] <- inverse(Sigma[,])
  }
  ", file="jagsmodel.txt", append=TRUE)
    
    #####################
    ## generate inits
    #####################
    if(is.null(svals)){
        gen.inits <- function(i, seed){
            .tmp <- list(cuts0 = sort(seq(from = -1.5, to = 1.5, length=ncol(ngk)-1) + rnorm(ncol(ngk)-1, sd=0.1)),
                         epsilon = cbind(rnorm(nrow(ngk),sd=0.2), rnorm(nrow(ngk),sd=0.1)),
                         eps_sd = runif(2,0,1),
                         eps_corr = runif(1,-1,1),
                         .RNG.name = "base::Mersenne-Twister",
                         .RNG.seed = seed)
            if(meanX){
                .tmp$beta_m <- rnorm(dat$dimXm, sd=0.0)
            }
            if(sdX){
                .tmp$beta_s <- rnorm(dat$dimXs, sd=0.0)
            }
            if(fixcuts){
                .tmp$mu0 <- rnorm(2, mean=0, sd=0.1)
                if(ncol(ngk) < 4){
                    .tmp$cuts0 <- NULL
                } else {
                    .tmp$cuts0 <- sort(seq(from = cuts12[2] + 1, to = cuts12[2] + 3, length=ncol(ngk)-3))
                }
            }
            dump(names(.tmp), envir=as.environment(.tmp), file=paste("jagsinits",i,".txt",sep=""))
        }
    } else { ## initialize using slightly jittered starting values
        gen.inits <- function(i, seed){
            .tmp <- list(cuts0     = sort(svals$cutpoints + rnorm(ncol(ngk)-1, sd=0.02)),
                         epsilon   = cbind(svals$mug, log(svals$sigmag)) + matrix(rnorm(2*nrow(ngk), sd=0.01), ncol=2),
                         eps_sd    = c(sd(svals$mug), sd(log(svals$sigmag))),
                         eps_corr  = cor(svals$mug, log(svals$sigmag)),
                         .RNG.name = "base::Mersenne-Twister",
                         .RNG.seed = seed)
            if(meanX){
                .tmp$beta_m      <- as.vector(coef(lm(svals$mug ~ dat$Xm - 1)))
                .tmp$epsilon[,1] <- as.vector(resid(lm(svals$mug ~ dat$Xm - 1)))
                .tmp$eps_sd[1]   <- sd(.tmp$epsilon[,1])
            }
            if(sdX){
                .tmp$beta_s      <- as.vector(coef(lm(I(log(svals$sigmag)) ~ dat$Xm - 1)))
                .tmp$epsilon[,2] <- as.vector(resid(lm(I(log(svals$sigmag)) ~ dat$Xm - 1)))
                .tmp$eps_sd[2]   <- sd(.tmp$epsilon[,2])        
            }
            dump(names(.tmp), envir=as.environment(.tmp), file=paste("jagsinits",i,".txt",sep=""))            
        }
    }
  
    gen.inits(1, bugs.seeds[1])
    gen.inits(2, bugs.seeds[2])
    gen.inits(3, bugs.seeds[3])

    ####################
    ## write script
    ####################
    cat("load dic\n", file="jagsscript.txt")
    if(useglm){
        cat("load glm\n", file="jagsscript.txt", append=TRUE)
    }
    cat("model in 'jagsmodel.txt'
data in 'jagsdata.txt'
compile, nchains(3)
parameters in 'jagsinits1.txt', chain(1)
parameters in 'jagsinits2.txt', chain(2)
parameters in 'jagsinits3.txt', chain(3)
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
    if(fixcuts){
        cat("monitor set mu0\n", file="jagsscript.txt", append=TRUE)
    }
    cat("monitor set eps_sd
monitor set eps_corr
monitor set pD
monitor set deviance\n", file="jagsscript.txt", append=TRUE)
    cat(paste("update ",niter,", by (20)\n",sep=""), file="jagsscript.txt", append=TRUE)
    cat("coda *, stem(jags)\n", file = "jagsscript.txt", append=TRUE)
  
    ############################
    ## run JAGS and get results
    ############################
    system("jags jagsscript.txt")
    ## kat change
    ##system('"C:/Program Files/JAGS/JAGS-4.0.0/x64/bin/jags.bat" jagsscript.txt') 
    ## system('"/home/research/jrlockwood/bin/jags" jagsscript.txt') 
  
    ## note: pD information is stored in its own CODA file because it pools across
    ## chains
    pD <- read.coda("jagschain0.txt","jagsindex0.txt", quiet=TRUE)
    r1 <- read.coda("jagschain1.txt", "jagsindex.txt", quiet=TRUE)
    r2 <- read.coda("jagschain2.txt", "jagsindex.txt", quiet=TRUE)
    r3 <- read.coda("jagschain3.txt", "jagsindex.txt", quiet=TRUE)
    rall <- rbind(r1, r2, r3)
    
    ## note deviance is at end, grab it and take it off and do DIC calc
    ## DIC = mean(deviance) + mean(pD)
    ## NOTE: from thread here:
    ## http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/ea46dc43/
    deviance <- rall[,"deviance"]
    rall <- rall[,setdiff(colnames(rall),"deviance")]
    dic <- mean(deviance) + mean(pD)
    
    gr <- vector(ncol(rall), mode="list")
    for(param in 1:ncol(rall)){
        gr[[param]] <- gelman.diag(mcmc.list(r1[,param],r2[,param], r3[,param]))$psrf
    }
    gr <- as.data.frame(do.call("rbind", gr))
    names(gr) <- c("gr.est","gr.upper")
    r <- data.frame(mean = apply(rall, 2, mean), sd = apply(rall, 2, sd), gr)
    cat("\n\nPARAMETERS WITH HIGH GR\n")
    print(subset(r, gr.est >= 1.05), digits=4)
    
    return(list(samps = rall, summaries = r, deviance = deviance, pD = pD, dic = dic))  
}

