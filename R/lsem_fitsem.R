## File Name: lsem_fitsem.R
## File Version: 0.5412

lsem_fitsem <- function( dat, weights, lavfit, fit_measures, NF, G, moderator.grid,
                verbose, pars, standardized, variables_model, sufficient_statistics,
                lavaan_fct, lavmodel, use_lavaan_survey=TRUE, pseudo_weights=0,
                est_joint=FALSE, par_invariant=NULL, par_linear=NULL, par_quadratic=NULL, 
                partable_joint=NULL, ... )
{

    parameters <- NULL
    fits <- NULL
    fitstats_joint <- NULL
    sample_stats <- NULL
    survey.fit <- lavfit
    pars0 <- pars
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)

    est_separate <- ! est_joint


    #- verbose
    pr <- lsem_fitsem_verbose_start(G=G, verbose=verbose)

    #- estimate model meanstructure
    is_meanstructure <- lavfit@Model@meanstructure

    #- sufficient statistics
    if (sufficient_statistics){
        sample_stats <- lsem_fitsem_compute_sufficient_statistics(G=G, dat=dat,
                    variables_model=variables_model, weights=weights)
    }

    #-- joint estimation
    if (est_joint){
        #- parameter table
        if ( is.null(partable_joint) ){
            partable_joint <- lsem_fitsem_joint_estimation_partable(lavfit=lavfit, G=G,
                                par_invariant=par_invariant, par_linear=par_linear,
                                par_quadratic=par_quadratic)
        }
        #- fit model
        if (sufficient_statistics){
            survey.fit <- lsem_fitsem_joint_estimation_sufficient_statistics(
                            partable_joint=partable_joint, is_meanstructure=is_meanstructure,
                            sample_stats=sample_stats, lavaan_est_fun=lavaan_est_fun, ... )
        }
        fitstats_joint <- lsem_lavaan_fit_measures(object=survey.fit, fit_measures=fit_measures)
    }

    #-- separate estimation: loop over groups
    for (gg in 1:G){
        dat$weight <- weights[,gg]

        #***** fit the model using weighted data
        if (( ! sufficient_statistics) & ( est_separate)){
            if (use_lavaan_survey){   # fit in lavaan.survey
                survey.fit <- lsem_fitsem_raw_data_lavaan_survey(dat=dat,
                                    lavmodel=lavmodel, lavfit=lavfit)
            }
            if (! use_lavaan_survey){  # fit in lavaan
                survey.fit <- lsem_fitsem_raw_data_lavaan(dat=dat, pseudo_weights=pseudo_weights,
                                    survey.fit=survey.fit, lavaan_est_fun=lavaan_est_fun, ...)
            }
        }
        #***** fit the model using sufficient statistics
        if (sufficient_statistics & est_separate){
            survey.fit <- lsem_fitsem_sufficient_statistics_lavaan( gg=gg, lavmodel=lavmodel,
                                lavaan_est_fun=lavaan_est_fun, survey.fit=survey.fit,
                                sample_stats=sample_stats, is_meanstructure=is_meanstructure, ... )
        }

        dfr.gg <- pars <- sirt_import_lavaan_parameterEstimates(object=survey.fit)

        if (standardized){
            sol <- sirt_import_lavaan_standardizedSolution(object=survey.fit)
            colnames(sol)[ which( colnames(sol)=="est.std" ) ] <- "est"
            sol$lhs <- paste0( "std__", sol$lhs)
            pars <- sirt_rbind_fill( x=pars, y=sol )
            dfr.gg <- pars
        }

        pars <- sirt_lavaan_partable_parnames(partable=pars)
        NP <- length(pars0)
        if (est_separate){
            ind <- match(pars0, pars)
            grid_index <- gg
            parindex <- 1:NP
            par_gg <- pars0
        } else {
            ind <- NULL
            # parindex <- NULL
            par_gg <- NULL
            for (pp in 1:NP){
                ind_pp <- which(pars==pars0[pp])
                npp <- length(ind_pp)
                ind <- c( ind, ind_pp)
                # parindex <- c(parindex, rep(pp, npp))
                par_gg <- c(par_gg, rep(pars0[pp], npp) )
            }
            grid_index <- dfr.gg$block[ ind ]
            parindex <- rep(1:NP, each=G)
            dfr.gg$block <- dfr.gg$group <- NULL
        }
        dfr.gg <- dfr.gg[ ind, ]
        dfr.gg <- data.frame(grid_index=grid_index, moderator=moderator.grid[grid_index],
                            par=par_gg, parindex=parindex, dfr.gg )
        if (est_separate){
            est_fit <- lsem_lavaan_fit_measures(object=survey.fit, fit_measures=fit_measures)
            dfr.gg0 <- data.frame(grid_index=gg, moderator=moderator.grid[gg],
                              par=fit_measures, parindex=NP+1:NF, est=est_fit, op='fit' )
            vars <- setdiff( colnames(dfr.gg), colnames(dfr.gg0) )
            for (vv in vars){ dfr.gg0[,vv] <- NA }
            dfr.gg <- rbind( dfr.gg, dfr.gg0[, colnames(dfr.gg) ] )
        }
        parameters <- rbind( parameters, dfr.gg )
        #- verbose
        res <- lsem_fitsem_verbose_progress(gg=gg, G=G, pr=pr, verbose=verbose)
        if (est_joint){ break }
    }

    res <- lsem_fitsem_verbose_progress(gg=G+1, G=G, pr=pr, verbose=verbose)
    parameters <- parameters[ order(parameters$parindex), ]
    rownames(parameters) <- paste0( parameters$par, "__", parameters$grid_index )

    #--- OUTPUT
    res <- list( parameters=parameters, is_meanstructure=is_meanstructure,
                    partable_joint=partable_joint, fitstats_joint=fitstats_joint,
                    sample_stats=sample_stats)
    return(res)
}

lsem.fitsem <- lsem_fitsem
