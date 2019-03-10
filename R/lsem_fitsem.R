## File Name: lsem_fitsem.R
## File Version: 0.478

lsem_fitsem <- function( dat, weights, lavfit,
            fit_measures, NF, G, moderator.grid, verbose,
            pars, standardized, variables_model,
            sufficient_statistics, lavaan_fct, lavmodel,
            use_lavaan_survey=TRUE, pseudo_weights=0, ... )
{

    parameters <- NULL
    fits <- NULL
    survey.fit <- lavfit
    pars0 <- pars
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)

    if (verbose){
        cat( "** Fit lavaan model\n")
        G1 <- min(G,10)
        pr <- round( seq(1,G, len=G1) )
        cat("|")
        cat( paste0( rep("*",G1), collapse="") )
        cat("|\n")
        cat("|")
    }

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
                # use starting values
                partable <- lavaan::parameterTable(object=survey.fit)
                partable$start <- partable$est
                survey.fit <- lavaan_est_fun(model=partable, data=dat1,
                                sampling.weights=sampling_weights, ... )
            }
        }
        #***** fit the model using sufficient statistics
        if (sufficient_statistics){
            res <- lsem_weighted_mean( x=dat[, variables_model], weights=dat$weight )
            wmean <- res$mean
            res <- lsem_weighted_cov( x=dat[, variables_model], weights=dat$weight )
            wcov <- res$cov
            Nobs <- round(res$Nobs)
            survey.fit <- lavaan_est_fun(model=lavmodel, sample.cov=wcov,
                                sample.mean=wmean, sample.nobs=Nobs, ... )
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
        dfr.gg <- data.frame("grid_index"=gg, "moderator"=moderator.grid[gg],
                          "par"=pars0, "parindex"=1:NP, dfr.gg    )
        est_fit <- lavaan::fitMeasures(object=survey.fit, fit.measures=fit_measures )
        dfr.gg0 <- data.frame("grid_index"=gg, "moderator"=moderator.grid[gg],
                          "par"=fit_measures, "parindex"=NP+1:NF,
                          "est"=est_fit, "op"="fit" )
        vars <- setdiff( colnames(dfr.gg), colnames(dfr.gg0) )
        for (vv in vars){ dfr.gg0[,vv] <- NA }
        dfr.gg <- rbind( dfr.gg, dfr.gg0[, colnames(dfr.gg) ] )
        parameters <- rbind( parameters, dfr.gg )
        # fits <- rbind( fits, dfr.gg )
        if (verbose){
            if ( gg %in% pr ){
                cat("-")
                utils::flush.console()
            }
        }
    }
    if (verbose){
        cat("|\n")
        utils::flush.console()
    }
    parameters <- parameters[ order(parameters$parindex), ]
    #--- OUTPUT
    res <- list( parameters=parameters )
    return(res)
}

lsem.fitsem <- lsem_fitsem
