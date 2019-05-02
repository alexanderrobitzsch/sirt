## File Name: lsem_fit_initial_model_sufficient_statistics.R
## File Version: 0.03

lsem_fit_initial_model_sufficient_statistics <- function(dat, variables_model,
    sampling_weights, has_meanstructure)
{
    data_suff <- dat[, variables_model]
    dat_resp <- 1 - is.na(data_suff)
    wmean <- lsem_weighted_mean( x=data_suff, weights=sampling_weights, x_resp=dat_resp)$mean
    res <- lsem_weighted_cov( x=data_suff, weights=sampling_weights, x_resp=dat_resp)
    wcov <- res$cov
    Nobs <- round(res$Nobs)
    if (! has_meanstructure){
        wmean <- NULL
    }
    #--- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs )
    return(res)
}
