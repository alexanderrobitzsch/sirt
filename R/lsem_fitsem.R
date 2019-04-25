## File Name: lsem_fitsem.R
## File Version: 0.502

lsem_fitsem <- function( dat, weights, lavfit, fit_measures, NF, G, moderator.grid,
                verbose, pars, standardized, variables_model, sufficient_statistics,
                lavaan_fct, lavmodel, use_lavaan_survey=TRUE, pseudo_weights=0, ... )
{

    parameters <- NULL
    fits <- NULL
    survey.fit <- lavfit
    pars0 <- pars
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)

    #- verbose
    pr <- lsem_fitsem_verbose_start(G=G, verbose=verbose)

    #- sufficient statistics
    if (sufficient_statistics){
        wmean <- wcov <- Nobs <- as.list(1:G)
        data_suff <- dat[, variables_model]
        dat_resp <- 1 - is.na(data_suff)
        for (gg in 1:G){
            weights_gg <- weights[,gg]
            res <- lsem_weighted_mean( x=data_suff, weights=weights_gg,    x_resp=dat_resp)
            wmean[[gg]] <- res$mean
            res <- lsem_weighted_cov( x=data_suff, weights=weights_gg, x_resp=dat_resp)
            wcov[[gg]] <- res$cov
            Nobs[[gg]] <- round(res$Nobs)
        }
    }

    #-- loop over groups
    for (gg in 1:G){
        dat$weight <- weights[,gg]

        #***** fit the model using weighted data
        if (! sufficient_statistics){
            if (use_lavaan_survey){
                survey.fit <- lsem_fitsem_raw_data_lavaan_survey(dat=dat,
                                    lavmodel=lavmodel, lavfit=lavfit)
            }
            if (! use_lavaan_survey){
                res <- lsem_fitsem_raw_data_define_pseudo_weights(dat=dat,
                            pseudo_weights=pseudo_weights)
                dat1 <- res$dat
                sampling_weights <- res$sampling_weights
                nobs_pseudo <- res$nobs_pseudo
                sum_weight <- res$sum_weight

                # use starting values
                partable <- lavaan::parameterTable(object=survey.fit)
                partable$start <- partable$est
                #- fit model
                survey.fit <- lavaan_est_fun(model=partable, data=dat1,
                                    sampling.weights=sampling_weights, ... )
                #- adjust sample size
                survey.fit <- lavaan_object_adjust_sample_size(object=survey.fit,
                                    n_used=sum_weight)
            }
        }
        #***** fit the model using sufficient statistics
        if (sufficient_statistics){
            survey.fit <- lavaan_est_fun(model=lavmodel, sample.cov=wcov[[gg]],
                                sample.mean=wmean[[gg]], sample.nobs=Nobs[[gg]], ... )
        }

        dfr.gg <- pars <- lavaan::parameterEstimates(object=survey.fit)
        if (standardized){
            sol <- lavaan::standardizedSolution(object=survey.fit)
            colnames(sol)[ which( colnames(sol)=="est.std" ) ] <- "est"
            sol$lhs <- paste0( "std__", sol$lhs)
            pars <- sirt_rbind_fill( x=pars, y=sol )
            dfr.gg <- pars
        }
        pars <- paste0( pars$lhs, pars$op, pars$rhs )
        NP <- length(pars0)
        ind <- match( pars0, pars )
        dfr.gg <- dfr.gg[ ind, ]
        dfr.gg <- data.frame(grid_index=gg, moderator=moderator.grid[gg],
                          par=pars0, parindex=1:NP, dfr.gg )
        est_fit <- lavaan::fitMeasures(object=survey.fit, fit.measures=fit_measures )
        dfr.gg0 <- data.frame(grid_index=gg, moderator=moderator.grid[gg],
                          par=fit_measures, parindex=NP+1:NF, est=est_fit, op='fit' )
        vars <- setdiff( colnames(dfr.gg), colnames(dfr.gg0) )
        for (vv in vars){ dfr.gg0[,vv] <- NA }
        dfr.gg <- rbind( dfr.gg, dfr.gg0[, colnames(dfr.gg) ] )
        parameters <- rbind( parameters, dfr.gg )
        #- verbose
        res <- lsem_fitsem_verbose_progress(gg=gg, G=G, pr=pr, verbose=verbose)
    }
    res <- lsem_fitsem_verbose_progress(gg=G+1, G=G, pr=pr, verbose=verbose)
    parameters <- parameters[ order(parameters$parindex), ]
    #--- OUTPUT
    res <- list( parameters=parameters )
    return(res)
}

lsem.fitsem <- lsem_fitsem
