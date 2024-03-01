## File Name: lsem_fitsem_compute_sufficient_statistics.R
## File Version: 0.124

lsem_fitsem_compute_sufficient_statistics <- function(G, dat, variables_model,
    weights, moderator_variable=NULL, loc_linear_smooth=NULL, moderator.grid=NULL,
    pd=FALSE, residualized_intercepts=NULL,    has_meanstructure=FALSE,
    residualize=TRUE, is_imputed=FALSE, Nimp=0, moderator=NULL)
{
    wmean <- wcov <- Nobs <- as.list(1:G)

    dat0 <- dat
    if (!is_imputed){
        dat0 <- list(dat0)
    }
    Nimp <- max(1, Nimp)
    Nobs_tt <- wcov_tt <- wmean_tt <- list()

    for (gg in 1L:G){

        for (ii in 1L:Nimp){
            if (!is_imputed){
                weights_gg <- weights[,gg]
            } else {
                weights_gg <- weights[[ii]][,gg, drop=TRUE]
            }
            dat <- dat0[[ii]]
            data_suff <- dat[, variables_model]
            dat_resp <- 1 - is.na(data_suff)
            if (is_imputed){
                moderator_variable <- dat[, moderator]
            }

            # res <- lsem_weighted_mean( x=data_suff, weights=weights_gg, x_resp=dat_resp)
            res <- lsem_weighted_cov( x=data_suff, weights=weights_gg, x_resp=dat_resp,
                        moderator_variable=moderator_variable,
                        loc_linear_smooth=loc_linear_smooth,
                        moderator_value=moderator.grid[gg], pd=pd,
                        residualized_intercepts=residualized_intercepts,
                        has_meanstructure=has_meanstructure, residualize=residualize)
            wmean_tt[[ii]] <- res$mean
            wcov_tt[[ii]] <- res$cov
            Nobs_tt[[ii]] <- res$Nobs
        }  # end ii

        wmean[[gg]] <- lsem_aggregate_statistics(x=wmean_tt)
        wcov[[gg]] <- lsem_aggregate_statistics(x=wcov_tt)
        Nobs[[gg]] <- round( lsem_aggregate_statistics(x=Nobs_tt) )

    }  # end gg

    #** adapt if mean structure is requested
    if ( has_meanstructure & residualize ){
        for (gg in 1L:G){
            wmean[[gg]] <- residualized_intercepts[gg,]
        }
    }

    #- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs)
    return(res)
}
