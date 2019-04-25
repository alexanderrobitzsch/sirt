## File Name: lsem.estimate.R
## File Version: 0.957

# estimate LSEM model
lsem.estimate <- function( data, moderator, moderator.grid,
        lavmodel, type="LSEM", h=1.1, residualize=TRUE,
        fit_measures=c("rmsea","cfi","tli","gfi","srmr"), standardized=FALSE,
        standardized_type="std.all", lavaan_fct="sem", sufficient_statistics=TRUE,
        use_lavaan_survey=FALSE, pseudo_weights=0, sampling_weights=NULL,
        eps=1E-8, verbose=TRUE, ... )
{
    CALL <- match.call()
    s1 <- Sys.time()

    lavaan.args <- list(...)
    if (standardized){
        if ( type=="MGM"){
            stop("standardized=TRUE cannot be applied for type='MGM'")
        }
    }

    #- process arguments
    res <- lsem_estimate_proc_args( lavaan.args=lavaan.args, sufficient_statistics=sufficient_statistics,
                pseudo_weights=pseudo_weights, lavmodel=lavmodel, data=data,
                use_lavaan_survey=use_lavaan_survey )
    sufficient_statistics <- res$sufficient_statistics
    use_lavaan_survey <- res$use_lavaan_survey
    variables_model <- res$variables_model
    use_pseudo_weights <- res$use_pseudo_weights
    variables_ordered <- res$variables_ordered

    # group moderator if type="MGM"
    out <- lsem_group_moderator( data=data, type=type, moderator.grid=moderator.grid,
                moderator=moderator, residualize=residualize, h=h )
    data <- out$data
    moderator.grouped <- out$moderator.grouped
    h <- out$h
    residualize <- out$residualize
    moderator.grid <- out$moderator.grid
    # residualize input data
    out <- lsem_residualize( data=data, moderator=moderator, moderator.grid=moderator.grid,
                lavmodel=lavmodel, h=h, residualize=residualize, eps=eps, verbose=verbose,
                sampling_weights=sampling_weights)
    G <- out$G
    data <- out$data
    weights <- out$weights
    residualized_intercepts <- out$residualized_intercepts
    N <- out$N
    bw <- out$bw
    moderator.density <- out$moderator.density
    sampling_weights <- out$sampling_weights
    no_sampling_weights <- out$no_sampling_weights
    data$index <- seq(1,N)

    # unweighted fit of lavaan model
    dat <- data
    lavmodel__ <- lavmodel

    #* fit lavaan model
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)
    lavfit <- lavaan_est_fun(model=lavmodel__, data=dat,  ... )

    # extract variables which are in model and data frame
    pars <- lavaan::parameterEstimates(object=lavfit)
    fM <- lavaan::fitMeasures(object=lavfit)
    fit_measures <- intersect( fit_measures, names(fM) )    
    fitstat <- fM[ fit_measures ]
    NF <- length(fit_measures)
    
    #-- testing purposes
    if (FALSE){
        partable <- lavaan::parameterTable(lavfit)
        partable$start <- partable$est
        partable$free <- 0
        lavfit2 <- lavaan_est_fun(partable, data=dat, ... )
        pars <- lavaan::parameterEstimates(object=lavfit2)
        lavfit2 <- lsem_lavaan_modify_lavaan_object_test(object=lavfit2)
        fitstat2 <- lsem_lavaan_fit_measures(object=lavfit2, fit_measures=fit_measures)
    }

    if (standardized){
        sol <- lavaan::standardizedSolution( object=lavfit, type=standardized_type)
        colnames(sol)[ which( colnames(sol)=="est.std" ) ] <- "est"
        sol$lhs <- paste0( "std__", sol$lhs)
        pars <- sirt_rbind_fill( x=pars, y=sol )
    }
    pars <- apply( pars[, c("lhs", "op", "rhs" ) ], 1, FUN=function(ll){
                        paste0( ll[1], ll[2], ll[3] ) } )
    # fit LSEM for all moderator groups
    out2 <- lsem_fitsem( dat=dat, weights=weights, lavfit=lavfit,
                    fit_measures=fit_measures, NF=NF, G=G, moderator.grid=moderator.grid,
                    verbose=verbose, pars=pars,    standardized=standardized,
                    variables_model=variables_model, sufficient_statistics=sufficient_statistics,
                    lavaan_fct=lavaan_fct, lavmodel=lavmodel, use_lavaan_survey=use_lavaan_survey,
                    pseudo_weights=pseudo_weights, ... )
    parameters <- out2$parameters
    rownames(parameters) <- paste0( parameters$par, "__", parameters$grid_index )

    #**** parameter and fit statistics summary
    parameters_summary <- lsem_parameter_summary( parameters=parameters,
                                moderator.density=out$moderator.density, verbose=verbose )
    out$moderator.density$Neff <- colSums(weights)

    obji0 <- obji <- out$moderator.density
    obji$moderator <- obji$moderator
    obji$wgt <- obji$wgt
    obji$Neff <- obji$Neff
    dfr <- data.frame( M=colMeans( obji0[,-1] ), SD=apply( obji0[,-1], 2, stats::sd ),
                        min=apply( obji0[,-1], 2, min ), max=apply( obji0[,-1], 2, max ))
    dfr0 <- data.frame(M=mean( data[,moderator], na.rm=TRUE ), SD=out$sd.moderator,
                        min=min( data[, moderator], na.rm=TRUE ),
                        max=max( data[, moderator], na.rm=TRUE ) )
    obji <- rbind( dfr0, dfr )
    rownames(obji) <- NULL
    moderator.stat <- data.frame(variable=c("moderator","wgt", "Neff"), obji )

    # output
    s2 <- Sys.time()
    res <- list( parameters=parameters, weights=weights,
                    parameters_summary=parameters_summary,
                    bw=out$bw, h=h, N=out$N, moderator.density=out$moderator.density,
                    moderator.stat=moderator.stat, moderator.grouped=moderator.grouped,
                    m.moderator=mean( data[,moderator], na.rm=TRUE ),
                    sd.moderator=out$sd.moderator, moderator=moderator,
                    moderator.grid=moderator.grid, lavmodel=lavmodel, residualize=residualize,
                    data=data, residualized.intercepts=residualized_intercepts,
                    lavaan.args=lavaan.args, fit_measures=fit_measures, s1=s1, s2=s2,
                    standardized=standardized, standardized_type=standardized_type,
                    lavaan_fct=lavaan_fct, use_lavaan_survey=use_lavaan_survey,
                    pseudo_weights=pseudo_weights, use_pseudo_weights=use_pseudo_weights,
                    sufficient_statistics=sufficient_statistics,
                    variables_ordered=variables_ordered, sampling_weights=sampling_weights,
                    no_sampling_weights=no_sampling_weights, type=type, CALL=CALL )
    class(res) <- "lsem"
    return(res)
}
