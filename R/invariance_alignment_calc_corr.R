## File Name: invariance_alignment_calc_corr.R
## File Version: 0.03


# auxiliary function for calculation of correlations
invariance_alignment_calc_corr <- function(parsM)
{
    cM <- stats::cor(parsM)
    I <- ncol(cM)
    rbar <- ( sum(cM) - I )/ ( I^2 - I)
    return(rbar)
}
