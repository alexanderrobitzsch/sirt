## File Name: ccov_np_compute_ccov_sum_score.R
## File Version: 0.10

ccov_np_compute_ccov_sum_score <- function(score, data, use_rcpp=TRUE)
{
    scores <- sort(unique(score))
    wgt_score <- sirt_sum_norm(table(score))
    NS <- length(scores)
    ccov_ff <- rep(NA,NS)
    if (!use_rcpp){
        for (ss in 1:NS){
            i1 <- which(score==scores[ss])
            s1 <- stats::cov.wt(x=data[i1,], method="ML")
            ccov_ff[ss] <- s1$cov[1,2]
        }
    } else {
        index <- match(score, scores)-1
        ccov_ff <- sirt_rcpp_ccov_np_compute_ccov_sum_score( index=index,
                            NS=NS, data=as.matrix(data) )
    }
    ccov_aggr <- sum(wgt_score*ccov_ff)
    res <- list(ccov_ff=ccov_ff, scores=scores, ccov_aggr=ccov_aggr)
    return(res)
}
