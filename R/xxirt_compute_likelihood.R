## File Name: xxirt_compute_likelihood.R
## File Version: 0.22


###########################################################################
xxirt_compute_likelihood <- function( probs_items, dat, resp_index,
        dat_resp_bool=NULL )
{
    N <- nrow(dat)
    TP <- dim(probs_items)[3]
    I <- dim(probs_items)[1]
    maxK <- dim(probs_items)[2]
    p.xi.aj <- matrix( 1, nrow=N, ncol=TP )
    probs <- matrix( probs_items, nrow=I, ncol=maxK*TP )
    p.xi.aj <- sirt_rcpp_xxirt_compute_likelihood( dat=dat, dat_resp_bool=dat_resp_bool,
                    probs=probs, TP=TP, maxK=maxK )
    return(p.xi.aj)
}
###########################################################################
