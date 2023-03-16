## File Name: lsdm_irf_distance_mad.R
## File Version: 0.04

lsdm_irf_distance_mad <- function(data, data_fitted, wgt_theta, use_abs=TRUE)
{
    dist <- data - data_fitted
    if (use_abs){
        dist <- abs(dist)
    }
    res <- rowSums(dist*wgt_theta)
    return(res)
}
