## File Name: lsem_fitsem_compute_sufficient_statistics.R
## File Version: 0.01

lsem_fitsem_compute_sufficient_statistics <- function(G, dat, variables_model,
    weights)
{
    wmean <- wcov <- Nobs <- as.list(1:G)
    data_suff <- dat[, variables_model]
    dat_resp <- 1 - is.na(data_suff)
    for (gg in 1:G){
        weights_gg <- weights[,gg]
        res <- lsem_weighted_mean( x=data_suff, weights=weights_gg, x_resp=dat_resp)
        wmean[[gg]] <- res$mean
        res <- lsem_weighted_cov( x=data_suff, weights=weights_gg, x_resp=dat_resp)
        wcov[[gg]] <- res$cov
        Nobs[[gg]] <- round(res$Nobs)
    }
    #- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs)
    return(res)
}
