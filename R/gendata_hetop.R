##########################################################################################
## function to generate count data by group from multinomial model with the
## structured probabilities under the HETOP
##########################################################################################
gendata_hetop <- function(G, K, Ng, mug, sigmag, cutpoints){
  stopifnot( (length(mug) == G) && length(sigmag)==G && length(Ng)==G)
  ## stopifnot(mug[1]==0 && sigmag[1]==1)
  stopifnot(length(cutpoints)==(K-1))
  
  ## calculate cell probabilities
  pgk <- matrix(0, ncol=K, nrow=G)
  for(g in 1:G){
    tmp     <- c(pnorm(cutpoints, mean=mug[g], sd = sigmag[g]), 1)
    pgk[g,] <- c(tmp[1], diff(tmp))
  }
  stopifnot(max(abs(apply(pgk, 1, sum) - 1)) < 1e-10)
  
  ## generate and return multinomial counts
  t(apply(cbind(Ng, pgk), 1,  function(x){ table(sample(as.factor(1:K), size=x[1], prob=x[-1], replace=T))}))
}
