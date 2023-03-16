## File Name: mcmcirt_compute_mh.R
## File Version: 0.04
## File Last Change: 2022-09-15


mcmcirt_compute_mh <- function(ll_old, ll_new)
{
    NL <- length(ll_new)
    mh <- ll_new - ll_old
    prob_mh <- ifelse(mh<0, exp(mh), 1)
    accept <- stats::runif(NL) < prob_mh
    #--- output
    ll_recent <- ifelse(accept, ll_new, ll_old)
    res <- list(accept=accept, ll_recent=ll_recent)
    return(res)
}
