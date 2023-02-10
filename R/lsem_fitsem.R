## File Name: lsem_fitsem.R
## File Version: 0.622

lsem_fitsem <- function( dat, weights, lavfit, fit_measures, NF, G, moderator.grid,
                verbose, pars, standardized, variables_model, sufficient_statistics,
                lavaan_fct, lavmodel, use_lavaan_survey=TRUE, pseudo_weights=0,
                est_joint=FALSE, par_invariant=NULL, par_linear=NULL, par_quadratic=NULL,
                partable_joint=NULL, pw_linear=1, pw_quadratic=1,
                se="standard", moderator_variable=NULL,
                loc_linear_smooth=NULL, pd=FALSE, has_meanstructure=FALSE,
                est_DIF=FALSE, ... )
{
    parameters <- NULL
    fits <- NULL
    fitstats_joint <- NULL
    sample_stats <- NULL
    data_joint <- NULL
    survey.fit <- lavfit
    pars0 <- pars
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)
    est_separate <- ! est_joint
    dif_effects <- NULL

    #- verbose
    pr <- lsem_fitsem_verbose_start(G=G, verbose=verbose)

    #- estimate model meanstructure
    is_meanstructure <- lavfit@Model@meanstructure

    #- sufficient statistics
    if (sufficient_statistics){
        sample_stats <- lsem_fitsem_compute_sufficient_statistics(G=G, dat=dat,
                    variables_model=variables_model, weights=weights,
                    moderator_variable=moderator_variable,
                    loc_linear_smooth=loc_linear_smooth, moderator.grid=moderator.grid,
                    pd=pd)
    }
    if (est_joint & (! sufficient_statistics)){
        N <- nrow(dat)
        ind <- rep(1:N, G)
        data_joint <- dat[ind, ]
        data_joint$group__ <- rep(1:G, each=N)
        data_joint$weight <- as.vector(weights)
        data_joint$case__ <- ind
    }

    #-- joint estimation
    if (est_joint){
        #- parameter table
        if ( is.null(partable_joint) ){
            partable_joint <- lsem_fitsem_joint_estimation_partable(lavfit=lavfit, G=G,
                                par_invariant=par_invariant, par_linear=par_linear,
                                par_quadratic=par_quadratic, pw_linear=pw_linear,
                                pw_quadratic=pw_quadratic)
        }
        #- fit model
        fit_args <- list( partable_joint=partable_joint,
                            is_meanstructure=is_meanstructure, data_joint=data_joint,
                            sample_stats=sample_stats, lavaan_est_fun=lavaan_est_fun,
                            verbose=verbose, se=se,
                            sufficient_statistics=sufficient_statistics,
                            G=G, pseudo_weights=pseudo_weights, ... )
        survey.fit <- do.call(what=lsem_fitsem_joint_estimation, args=fit_args)
        partable_joint <- lavaan::parameterTable(object=survey.fit)
        fitstats_joint <- lsem_lavaan_fit_measures(object=survey.fit,
                                fit_measures=fit_measures)

        #*** assessment of DIF effects
        partable_dif <- partable_joint
        partable_dif[ partable_dif$free > 0, "free"] <- 0
        partable_dif[ partable_dif$par %in% par_invariant, "free"] <- 1
        ind1 <- intersect( grep("=~", partable_dif$par), which( partable_joint$free==0) )
        ind2 <- NULL
        if (has_meanstructure){
            ind2 <- intersect( grep("~1", partable_dif$par),
                                    which( partable_joint$free==0) )
        }
        ind3 <- intersect( grep("~~", partable_dif$par), which( partable_joint$free==0) )
        partable_dif[ c(ind1,ind2, ind3), "free"] <- 1

        partable_dif <- partable_dif[ - grep( "_con", partable_dif$par), ]
        partable_dif$start <- partable_dif$est
        partable_dif$free <- cumsum(partable_dif$free)*(partable_dif$free>0)
        if ( sum(partable_dif$free)==0 ){
            est_DIF <- FALSE
        }
        if ( est_DIF ){
            fit_args$partable_joint <- partable_dif
            fit_dif <- do.call(what=lsem_fitsem_joint_estimation, args=fit_args)
            dif_effects <- lavaan::parameterTable(object=fit_dif)
            dif_effects <- dif_effects[ order(dif_effects$par), ]
            dif_effects$DIF <- 1*( dif_effects$free > 0 )
        }

        #--- output joint estimation
        fitstats_joint <- data.frame(stat=names(fitstats_joint), value=fitstats_joint)
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
                survey.fit <- lsem_fitsem_raw_data_lavaan(dat=dat,
                                    pseudo_weights=pseudo_weights,
                                    survey.fit=survey.fit, lavaan_est_fun=lavaan_est_fun,
                                    se=se, ...)
            }
        }
        #***** fit the model using sufficient statistics
        if (sufficient_statistics & est_separate){
            survey.fit <- lsem_fitsem_sufficient_statistics_lavaan( gg=gg,
                                lavmodel=lavmodel, lavaan_est_fun=lavaan_est_fun,
                                survey.fit=survey.fit, sample_stats=sample_stats,
                                is_meanstructure=is_meanstructure, se=se, ... )
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
            par_gg <- NULL
            for (pp in 1:NP){
                ind_pp <- which(pars==pars0[pp])
                npp <- length(ind_pp)
                ind <- c( ind, ind_pp)
                par_gg <- c(par_gg, rep(pars0[pp], npp) )
            }
            grid_index <- dfr.gg$group[ ind ]
            parindex <- rep(1:NP, each=G)
            dfr.gg$block <- dfr.gg$group <- NULL
        }
        dfr.gg <- dfr.gg[ ind, ]
        dfr.gg <- data.frame(grid_index=grid_index, moderator=moderator.grid[grid_index],
                            par=par_gg, parindex=parindex, dfr.gg )
        if (est_separate){
            est_fit <- lsem_lavaan_fit_measures(object=survey.fit,
                                fit_measures=fit_measures)
            dfr.gg0 <- data.frame(grid_index=gg, moderator=moderator.grid[gg],
                              par=fit_measures, parindex=NP+1:NF, est=est_fit, op='fit' )
            vars <- setdiff( colnames(dfr.gg), colnames(dfr.gg0) )
            for (vv in vars){ dfr.gg0[,vv] <- NA }
            dfr.gg <- rbind( dfr.gg, dfr.gg0[, colnames(dfr.gg) ] )
        }
        parameters <- rbind( parameters, dfr.gg )

        if (est_DIF){
            dif1 <- dif_effects[ dif_effects$DIF==1, ]
            pars1 <- intersect(paste(dif1$par), paste(parameters$par))
            dif1 <- dif1[ dif1$par %in% pars1, ]
            NP <- max(parameters$parindex)
            parameters1 <- parameters[ parameters$par %in% pars1, ]
            parameters1$est <- dif1$est
            parameters1$par <- paste0("dif__", parameters1$par)
            parameters1$parindex <- match( parameters1$parindex,
                                            unique( parameters1$parindex ) ) + NP
            parameters <- rbind( parameters, parameters1)
        }


        #- verbose
        res <- lsem_fitsem_verbose_progress(gg=gg, G=G, pr=pr, verbose=verbose)
        if (est_joint){ break }
    }
    if ( sum(colnames(parameters)=="se")==0){
        parameters$se <- NA
    }

    res <- lsem_fitsem_verbose_progress(gg=G+1, G=G, pr=pr, verbose=verbose)
    parameters <- parameters[ order(parameters$parindex), ]
    rownames(parameters) <- paste0( parameters$par, "__", parameters$grid_index )

    #--- OUTPUT
    res <- list( parameters=parameters, is_meanstructure=is_meanstructure,
                    partable_joint=partable_joint, fitstats_joint=fitstats_joint,
                    sample_stats=sample_stats, dif_effects=dif_effects)
    return(res)
}

lsem.fitsem <- lsem_fitsem
