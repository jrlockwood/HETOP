##########################################################################################
## function to find a cutpoint corresponding to a given cumulative probability
##########################################################################################
findcut <- function(mug, sigmag, pg, targetprob){
    stopifnot( (length(mug) == length(sigmag)) && (length(mug) == length(pg)) )
    stopifnot( max(abs(sum(pg) - 1)) < 1e-7)
    stopifnot( (targetprob > 0) & (targetprob < 1) )

    f <- function(theta){
        (sum(pg * pnorm(theta, mean = mug, sd = sigmag)) - targetprob)^2
    }
    return(nlm(f, 0.0)$estimate)
}
