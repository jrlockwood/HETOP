##########################################################################################
## function to compute grids from {ngk} and two fixed cuts
##########################################################################################
compute_grids <- function(ngk, fixedcuts){
    stopifnot(fixedcuts[1] < fixedcuts[2])
    G <- nrow(ngk)
    K <- ncol(ngk)

    ## collapse data to three categories
    ngkcollapse <- cbind(ngk[,1], ngk[,2], apply(ngk[,3:K,drop=F], 1, sum))

    ## if group has any zero cells, add 0.5 to all cells for that group
    zeros <- which(apply(ngkcollapse, 1, function(x){ any(x==0) }))
    ngkcollapse[zeros,] <- ngkcollapse[zeros,] + 0.5

    ## get MLE of (mu, logsigma) conditional on fixedcuts for each group, along with asymptotic SEs
    vals <- list()
    for(g in 1:nrow(ngkcollapse)){
        negll <- function(param){
            tmp    <- c(pnorm(fixedcuts, mean= param[1], sd = exp(param[2])), 1)
            pgk    <- c(tmp[1], diff(tmp))
            pgk[which(pgk <= 1e-300)]     <- 1e-300
            pgk[which(pgk >= 1 - 1e-300)] <- 1 - 1e-300
            -sum(ngkcollapse[g,] * log(pgk))
        }
        mustart <- weighted.mean(c(-1,0,1), w = ngkcollapse[g,])
        res <- nlm(f = negll, hessian=TRUE, print.level=0, p = c(mustart,0), iterlim=1500)
        acov <- solve(res$hessian)
        
        vals[[g]] <- data.frame(muhat    = res$estimate[1],
                                lshat    = res$estimate[2],
                                muhat_se = sqrt(acov[1,1]),
                                lshat_se = sqrt(acov[2,2]))
    }
    vals <- do.call("rbind", vals)

    res <- list()
    
    ## get extreme estimates
    res$mu.hatrange <- range(vals$muhat)
    res$ls.hatrange <- range(vals$lshat)

    ## get ranges using mean +/- 4*estimated std dev.
    ## estimate std dev using MOM and trimming to avoid ridiculous values
    mu.sd           <- sqrt(var(vals$muhat) - mean(vals$muhat_se^2, trim=0.01))
    res$mu.momrange <- c(mean(vals$muhat) - 4*mu.sd, mean(vals$muhat) + 4*mu.sd)

    ls.sd           <- sqrt(var(vals$lshat) - mean(vals$lshat_se^2, trim=0.01))
    res$ls.momrange <- c(mean(vals$lshat) - 4*ls.sd, mean(vals$lshat) + 4*ls.sd)
    return(res)
}
