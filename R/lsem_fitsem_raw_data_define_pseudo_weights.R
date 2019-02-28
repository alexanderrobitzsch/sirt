## File Name: lsem_fitsem_raw_data_define_pseudo_weights.R
## File Version: 0.01

lsem_fitsem_raw_data_define_pseudo_weights <- function(dat, pseudo_weights)
{
    sampling_weights <- "weight"
    if (pseudo_weights>0){
        weights <- dat$weight
        weights <- round( weights * pseudo_weights )
        N <- nrow(dat)
        ind <- rep(1:N, weights)
        dat <- dat[ind,]
        sampling_weights <- NULL
    }
    #-- output
    res <- list(dat=dat, sampling_weights=sampling_weights)
    return(res)
}
