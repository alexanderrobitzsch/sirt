## File Name: rm_proc_create_pseudoraters.R
## File Version: 0.10

rm_proc_create_pseudoraters <- function( dat, rater, pid, reference_rater=NULL )
{
    dat0 <- dat
    rater0 <- paste0(rater)
    pid0 <- pid
    items <- colnames(dat)
    dat <- NULL
    pid <- NULL
    rater <- NULL
    I <- length(items)

    m0 <- as.data.frame( matrix(NA, nrow=nrow(dat0), ncol=I) )
    colnames(m0) <- colnames(dat0)

    for (ii in 1:I){
        dat_ii <- m0
        dat_ii[, ii ] <- dat0[,ii]
        rater_ii <- paste0( rater0, "-", items[ii] )
        rater <- c( paste(rater), paste(rater_ii) )
        dat <- rbind( dat, dat_ii)
        pid <- c( pid, pid0)
    }

    if ( ! is.null(reference_rater) ){
        reference_rater <- paste0(reference_rater, "-", items )
    }

    #-- output
    res <- list( dat=dat, rater=rater, pid=pid, reference_rater=reference_rater)
    return(res)
}
