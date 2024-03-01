## File Name: lsem_fit_initial_model_sufficient_statistics.R
## File Version: 0.075

lsem_fit_initial_model_sufficient_statistics <- function(dat, variables_model,
    sampling_weights, has_meanstructure, is_imputed=FALSE, Nimp=0)
{

    dat0 <- dat
    if (!is_imputed){
        dat0 <- list(dat0)
    }
    Nimp <- max(1, Nimp)
    wmean <- list()
    wcov <- list()
    Nobs <- list()
    for (ii in 1L:Nimp){
        dat <- dat0[[ii]]
        data_suff <- dat[, variables_model]
        dat_resp <- 1 - is.na(data_suff)
        wmean[[ii]] <- lsem_weighted_mean( x=data_suff, weights=sampling_weights,
                                        x_resp=dat_resp)$mean
        res <- lsem_weighted_cov( x=data_suff, weights=sampling_weights, x_resp=dat_resp)
        wcov[[ii]] <- res$cov
        Nobs[[ii]] <- round(res$Nobs)
    }

    Nobs <- lsem_aggregate_statistics(x=Nobs)
    wcov <- lsem_aggregate_statistics(x=wcov)
    wmean <- lsem_aggregate_statistics(x=wmean)

    if (! has_meanstructure){
        wmean <- NULL
    }

    #--- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs )
    return(res)
}
