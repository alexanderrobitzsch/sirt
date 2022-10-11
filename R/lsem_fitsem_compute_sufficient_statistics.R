## File Name: lsem_fitsem_compute_sufficient_statistics.R
## File Version: 0.094

lsem_fitsem_compute_sufficient_statistics <- function(G, dat, variables_model,
    weights, moderator_variable=NULL, loc_linear_smooth=NULL, moderator.grid=NULL,
    pd=FALSE)
{
    wmean <- wcov <- Nobs <- as.list(1:G)
    data_suff <- dat[, variables_model]
    dat_resp <- 1 - is.na(data_suff)
    for (gg in 1:G){
        weights_gg <- weights[,gg]
        # res <- lsem_weighted_mean( x=data_suff, weights=weights_gg, x_resp=dat_resp)
        res <- lsem_weighted_cov( x=data_suff, weights=weights_gg, x_resp=dat_resp,
                    moderator_variable=moderator_variable,
                    loc_linear_smooth=loc_linear_smooth,
                    moderator_value=moderator.grid[gg], pd=pd)
        wmean[[gg]] <- res$mean
        wcov[[gg]] <- res$cov
        Nobs[[gg]] <- round(res$Nobs)
    }

    #- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs)
    return(res)
}
