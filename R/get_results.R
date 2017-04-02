get_results <- function(rep, data) {
  # data: matrix of counts, 1 row per group, column per category; REPS stacked
  # rep: which rows of data to use (which REP)
  # changes for 3/4/16:
  #   try to estimate MLEs even if "bad case"
  #   compare the bias-corrected and actual MLEs
  #   use the new est_deref_jags() function
  
  dat <- as.matrix(data[ data$rep == rep , c('level0', 'level1', 'level2', 'level3')])
  colnames(dat) <- c('1', '2', '3', '4')
  rownames(dat) <- NULL
  Ng <- apply(dat, 1, sum)
  Pg <- Ng/sum(Ng)
  G <- nrow(dat)
  K <- ncol(dat)
  
  ests <- list()
  
  # return NA for MLE estimates if a group has all observations in 1 cell
  # let's see what values we get if we take out this check
  #   chk <- apply(dat, 1, function(x) {sum(x==0)})
  #   if (sum(chk>2)>0) {
  #     ests <- data.frame(mlePrime = rep(NA, times = 2*G+1),
  #                        mleStar = rep(NA, times = 2*G+1))
  #   } else {
  #     # MLEs
  #     mle <- getmle_stz(dat, acov=FALSE)
  #     ests <- data.frame(mlePrime = c(NA, unlist(mle[c("mug","sigmag")])),
  #                        mleStar = unlist(est_deref(mprime = mle$mug, sprime = mle$sigmag, ngs = Ng, geticc = TRUE)))
  #   }
  #   rm(chk)
  
  mle <- getmle_stz(dat, acov=FALSE)
  
  mleicc <- sum( Pg * (mle$mug - sum(Pg * mle$mug))^2 ) / ( sum( Pg * (mle$mug - sum(Pg * mle$mug))^2 ) + sum( Pg * (mle$sigmag^2) ) )
  ests$mlePrime <- c(mleicc, unlist(mle[c("mug", "sigmag", "cutpoints")]))
  names(ests$mlePrime) <- c('icc', paste0("mug", 1:G), paste0("sigmag", 1:G), paste0("cutpoints", 1:(K-1)))
  
  ests$mleStar_biascor <- c(unlist(est_deref(mprime = mle$mug, sprime = mle$sigmag, cutests = mle$cutpoints, ngs = Ng, geticc = T, biascor=TRUE)))
  ests$mleStar <- c(unlist(est_deref(mprime = mle$mug, sprime = mle$sigmag, cutests = mle$cutpoints, ngs = Ng, geticc = T, biascor=FALSE)))
  names(ests$mleStar_biascor) <- names(ests$mleStar) <- c('icc', paste0("mstar", 1:G), paste0("sstar", 1:G), paste0("cutstar", 1:(K-1)))
  
  ests$mle_nlmecode <- mle$nlmdetails$code
  ests$mle_nlmeiter <- mle$nlmdetails$iterations
  
  # F-H's
  ok <- FALSE
  while(!ok){
    print(seedpos <- seedpos + 1)
    jags <- hetop_s1_fh_jags_v2f(dat, bugs.seeds=bugs.seeds[3*(seedpos-1)+c(1,2,3)], nburn=10000, niter=15000)
    if(nrow(subset(jags$summaries, gr.est > 1.1)) < 5){
      ok <- TRUE
    }
  }
  
  r <- jags$summaries
  print(r[c("cuts[1]","cuts[2]","cuts[3]","eps_sd[1]","eps_sd[2]","eps_corr"), ])
  jagsicc <- mean( apply(jags$samps[ ,1:(2*G)], 1, function(x) {
    m <- x[1:G]
    s <- x[(G+1):(2*G)]
    return( sum(Pg * (m - sum(Pg*m))^2) / (sum(Pg * (m - sum(Pg*m))^2) + sum(Pg * s^2)))
  }) )
  ests$jagsPrime <- c(jagsicc, r[c(grep("mu",rownames(r)), grep("sigma",rownames(r)), grep("cuts", rownames(r))),"mean"])
  names(ests$jagsPrime) <- c('icc', paste0("mug", 1:G), paste0("sigmag", 1:G), paste0("cutpoints", 1:(K-1)))
  
  # compute standardized "star" estimates at every REP of the MCMC chain
  ests$jagsStar <- apply( apply(jags$samps[,1:(2*G+(K-1))], 1,
                                function(X){ 
                                  unlist(est_deref_jags(mprime = X[(1:G)], sprime = X[(G+1):(2*G)], cutests = X[(2*G+1):(2*G+(K-1))], ngs = Ng, geticc = TRUE))
                                }),
                          1, mean)
  names(ests$jagsStar) <- c('icc', paste0("mstar", 1:G), paste0("sstar", 1:G), paste0("cutstar", 1:(K-1)))
  
  ests$jagsOther <- r[c("cuts[1]","cuts[2]","cuts[3]","eps_sd[1]","eps_sd[2]","eps_corr"), ]
  
  return(ests)
  
}
