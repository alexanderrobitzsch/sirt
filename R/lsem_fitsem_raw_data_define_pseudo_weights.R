## File Name: lsem_fitsem_raw_data_define_pseudo_weights.R
## File Version: 0.141

lsem_fitsem_raw_data_define_pseudo_weights <- function(dat, pseudo_weights)
{
    sampling_weights <- 'weight'
    W <- sum(dat[,sampling_weights])
    if (pseudo_weights>0){
        weights <- dat$weight
        N <- nrow(dat)
        W <- sum(weights)
        fac <- pseudo_weights / W
        weights1 <- fac*weights
        weights2 <- cumsum(weights1)
        tweights2 <- floor(weights2)
        freq <- c(tweights2[1], diff(tweights2))
        ind <- rep(1:N, freq)
        dat <- dat[ind,]
        sampling_weights <- NULL
    }
    nobs_pseudo <- nrow(dat)
    sum_weight <- W
    #-- output
    res <- list(dat=dat, sampling_weights=sampling_weights,
                    nobs_pseudo=nobs_pseudo, sum_weight=sum_weight)
    return(res)
}
