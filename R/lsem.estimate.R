## File Name: lsem.estimate.R
## File Version: 0.971

# estimate LSEM model
lsem.estimate <- function( data, moderator, moderator.grid,
        lavmodel, type="LSEM", h=1.1, bw=NULL, residualize=TRUE,
        fit_measures=c("rmsea","cfi","tli","gfi","srmr"), standardized=FALSE,
        standardized_type="std.all", lavaan_fct="sem", sufficient_statistics=FALSE,
        use_lavaan_survey=FALSE, pseudo_weights=0, sampling_weights=NULL,
        est_joint=FALSE, par_invariant=NULL, par_linear=NULL, par_quadratic=NULL,
        partable_joint=NULL, se=NULL, kernel="gaussian", eps=1E-8, verbose=TRUE, ... )
{
    lsem_args <- c(as.list(environment()), list(...))
    CALL <- match.call()
    s1 <- Sys.time()
    lavaan.args <- list(...)
    if (standardized){
        if ( type=="MGM"){
            stop("standardized=TRUE cannot be applied for type='MGM'")
        }
    }

    #- data cleaning
    data <- as.data.frame(data)
    data <- data[ ! is.na(data[,moderator]), ]

    #- process arguments
    res <- lsem_estimate_proc_args( lavaan.args=lavaan.args, sufficient_statistics=sufficient_statistics,
                pseudo_weights=pseudo_weights, lavmodel=lavmodel, data=data,
                use_lavaan_survey=use_lavaan_survey, est_joint=est_joint,
                par_invariant=par_invariant, par_linear=par_linear,
                par_quadratic=par_quadratic, partable_joint=partable_joint,
                moderator.grid=moderator.grid, se=se, verbose=verbose )
    sufficient_statistics <- res$sufficient_statistics
    use_lavaan_survey <- res$use_lavaan_survey
    variables_model <- res$variables_model
    use_pseudo_weights <- res$use_pseudo_weights
    variables_ordered <- res$variables_ordered
    est_joint <- res$est_joint
    partable <- res$partable
    has_meanstructure <- res$has_meanstructure
    se <- res$se
    compute_se <- res$compute_se
    pseudo_weights <- res$pseudo_weights

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
                lavmodel=lavmodel, h=h, bw=bw, residualize=residualize, eps=eps, verbose=verbose,
                sampling_weights=sampling_weights, kernel=kernel, variables_model=variables_model)
    G <- out$G
    data <- out$data
    weights <- out$weights
    residualized_intercepts <- out$residualized_intercepts
    N <- out$N
    bw <- out$bw
    h <- out$h
    moderator.density <- out$moderator.density
    sampling_weights <- out$sampling_weights
    no_sampling_weights <- out$no_sampling_weights
    m.moderator <- out$m.moderator
    sd.moderator <- out$sd.moderator
    data$index <- seq(1,N)

    # unweighted fit of lavaan model
    dat <- data
    lavmodel__ <- lavmodel

    #* extract estimation function
    lavaan_est_fun <- lsem_define_lavaan_est_fun(lavaan_fct=lavaan_fct)

    #* fit initial lavaan model
    lavfit <- lsem_fit_initial_model( lavmodel__=lavmodel__,
                    lavaan_est_fun=lavaan_est_fun, dat=dat, variables_model=variables_model,
                    sampling_weights=sampling_weights, has_meanstructure=has_meanstructure,
                    sufficient_statistics=sufficient_statistics, est_joint=est_joint,
                    se=se, use_lavaan_survey=use_lavaan_survey, ... )

    # extract variables which are in model and data frame
    pars <- sirt_import_lavaan_parameterEstimates(object=lavfit)
    fM <- sirt_import_lavaan_fitMeasures(object=lavfit)
    fit_measures <- intersect( fit_measures, names(fM) )
    fitstat <- fM[ fit_measures ]
    NF <- length(fit_measures)

    if (standardized){
        sol <- sirt_import_lavaan_standardizedSolution( object=lavfit, type=standardized_type)
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
                    pseudo_weights=pseudo_weights, est_joint=est_joint,
                    par_invariant=par_invariant, par_linear=par_linear,
                    par_quadratic=par_quadratic, partable_joint=partable_joint,
                    se=se, ... )
    parameters <- out2$parameters
    is_meanstructure <- out2$is_meanstructure
    fitstats_joint <- out2$fitstats_joint
    partable_joint <- out2$partable_joint
    sample_stats <- out2$sample_stats

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

    #-- model parameters
    model_parameters <- setdiff( paste(parameters_summary$par), fit_measures)

    #-- output
    s2 <- Sys.time()
    time <- s2-s1
    res <- list( parameters=parameters, weights=weights,
                    parameters_summary=parameters_summary,
                    bw=out$bw, h=h, N=out$N, moderator.density=out$moderator.density,
                    moderator.stat=moderator.stat, moderator.grouped=moderator.grouped,
                    m.moderator=m.moderator, sd.moderator=sd.moderator, moderator=moderator,
                    moderator.grid=moderator.grid, lavmodel=lavmodel, residualize=residualize,
                    data=data, residualized.intercepts=residualized_intercepts,
                    lavaan.args=lavaan.args, lsem_args=lsem_args, fit_measures=fit_measures,
                    model_parameters=model_parameters, s1=s1, s2=s2, time=time,
                    standardized=standardized, standardized_type=standardized_type,
                    lavaan_fct=lavaan_fct, use_lavaan_survey=use_lavaan_survey,
                    pseudo_weights=pseudo_weights, use_pseudo_weights=use_pseudo_weights,
                    sufficient_statistics=sufficient_statistics,
                    variables_ordered=variables_ordered, sampling_weights=sampling_weights,
                    no_sampling_weights=no_sampling_weights, is_meanstructure=is_meanstructure,
                    par_invariant=par_invariant, par_linear=par_linear, par_quadratic=par_quadratic,
                    est_joint=est_joint, fitstats_joint=fitstats_joint, partable_joint=partable_joint,
                    sample_stats=sample_stats, se=se, compute_se=compute_se,
                    class_boot=FALSE, type=type, CALL=CALL )
    class(res) <- "lsem"
    return(res)
}
