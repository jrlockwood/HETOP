###########################################################
## FUNCTION TO STANDARDIZE AND APPLY TRIPLE GOAL 
###########################################################

process_jags	<- function(ngk, Xm=NULL,Xs=NULL, bugs.seeds, nburn=5000, niter=3333,geticc=FALSE) {
  ##STEP 1: FH MODEL
  j1	<- hetop_s1_fh_jags_v2f(ngk=d, Xm=Xm,Xs=Xs, bugs.seeds, nburn=nburn, niter=niter)
  ##STEP 2: STANDARDIZE 
  G	<- nrow(ngk)
  tmp	<-  t(apply(j1$samps, 1,function(X){ 
    unlist(est_deref_jags(mprime = X[(1:G)], sprime = X[(G+1):(2*G)], cutests=X[(2*G+1):(2*G+3)], ngs = apply(ngk,1,sum), geticc = geticc))
  }))	
  ##STEP 3: APPLY TRIPLE GOAL TO MU*
  mu3	<- triple_goal(tmp[,grep("mstar",colnames(tmp))])
  ##STEP 4: APPLY TRIPLE GOAL TO SIGMA*
  sig3	<- triple_goal(tmp[,grep("sstar",colnames(tmp))])
  ##OUTPUT
  return(list(jagsobj=j1,starsamps=tmp,triple.mus=mu3,triple.sigs=sig3))
}
