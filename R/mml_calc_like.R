## File Name: mml_calc_like.R
## File Version: 0.04
## File Last Change: 2021-09-22


#-- calculation of the likelihood
mml_calc_like <- function (dat2, dat2resp, probs, pseudoll=0)
{
    if ( pseudoll==0 ){
        res <- MML2_CALCPOST_V1( DAT2=dat2, DAT2RESP=dat2resp, PROBS=probs)
    }
    if ( pseudoll==1 ){
        res <- sirt_rcpp_rasch_mml2_calcpost_pseudoll( DAT2=dat2, DAT2RESP=dat2resp,
                    PROBS=probs)
    }
    return(res)
}

