## File Name: lsem.permutationTest.R
## File Version: 0.593


#*** permutation test for LSEM model
lsem.permutationTest <- function( lsem.object, B=1000, residualize=TRUE,
            verbose=TRUE, n.core=1, cl.type="PSOCK" )
{

    s1 <- Sys.time()
    CALL <- match.call()

    lavaan.args <- lsem.object$lavaan.args
    entr <- c( 'lavmodel', 'data', 'h', 'moderator.grid', 'moderator', 'eps',
                    'fit_measures', 'par_invariant', 'est_joint',
                    'par_linear', 'par_quadratic', 'pw_linear',
                    'pw_quadratic')
    object <- lsem.object
    arglist <- list()
    EE <- length(entr)
    for (ee in 1L:EE){
        arglist[[ entr[ee] ]] <- object[[ entr[ee] ]]
    }
    arglist2 <- lsem.object$lavaan.args
    NL <- length(arglist2)
    if (NL > 0){
        for (ll in 1L:NL){
            arglist[[ names(arglist2)[ll] ]] <- arglist2[[ names(arglist2)[ll] ]]
        }
    }
    arglist$residualize <- residualize
    arglist$verbose <- FALSE
    arglist$standardized <- lsem.object$standardized
    arglist$standardized_type <- lsem.object$standardized_type
    arglist$lavaan_fct <- lsem.object$lavaan_fct
    arglist$sufficient_statistics <- lsem.object$sufficient_statistics
    use_lavaan_survey <- FALSE
    arglist$pseudo_weights <- lsem.object$pseudo_weights

    data1 <- data0 <- object$data
    parameters <- object$parameters
    parameters_summary <- object$parameters_summary
    moderator <- object$moderator

    #******************************************
    # start permutation test

    parameters_permutation <- matrix( NA, nrow(parameters), ncol=B)
    parameters_summary_M <- matrix( NA, nrow(parameters_summary), ncol=B)
    rownames(parameters_summary_M) <- parameters_summary$par
    parameters_summary_MAD <- parameters_summary_SD <- parameters_summary_M
    parameters_summary_lin_slo <- parameters_summary_M

    if ( verbose ){
        cat('Permutation test LSEM \n')
    }
    #** estimate model for each permutation
    bb <- 1
    nn <- 0

    N <- nrow(data0)

    # sampled values
    l1 <- list()

    #--- nonparallel version
    if (n.core==1){
        while (bb<=B){
            data1[, moderator ] <- sample( data0[, moderator ] )
            arglist$data <- data1
            res0 <- try( do.call( what='lsem.estimate', args=arglist ), silent=TRUE )
            res0_out <- list(est=res0$parameters$est, M=res0$parameters_summary$M,
                                SD=res0$parameters_summary$SD,
                                MAD=res0$parameters_summary$MAD,
                                lin_slo=res0$parameters_summary$lin_slo)
            res0_out <- lsem_permutationTest_collect_output(res0=res0)
            if ( ! inherits(res0,'try-error') ){
                parameters_permutation[, bb] <- res0_out$est
                parameters_summary_M[,bb] <- res0_out$M
                parameters_summary_SD[,bb] <- res0_out$SD
                parameters_summary_MAD[,bb] <- res0_out$MAD
                parameters_summary_lin_slo[,bb] <- res0_out$lin_slo
                if (verbose){
                    cat( bb, ' ')
                    if ( bb %% 20==0 ){
                        cat('\n')
                    }
                    utils::flush.console();
                }
                bb <- bb+1
            } else {
                nn <- nn + 1
                if (verbose){
                    cat( '\n Delete permutation dataset due to non-convergence\n')
                    utils::flush.console()
                }
            }
        }  # end bb
        if (verbose){ cat('\n') }
    }

    #------
    #* parallel version
    if (n.core>1){
        m <- parallel::detectCores()
        if(n.core>m-1){
            n.core <- m-1
        }
        cl <- parallel::makeCluster(spec=n.core, type=cl.type)
        varlist <- c('arglist')
        parallel::clusterExport(cl=cl, varlist=varlist, envir=environment() )
        parallel::clusterEvalQ(cl, expr=library(sirt))

        fun_lsem_permutation_test <- function(i, arglist)
        {
            data1 <- data0 <- arglist$data
            moderator <- arglist$moderator
            iterate <- TRUE
            while(iterate){
                data1[, moderator ] <- sample( data0[, moderator ] )
                arglist1 <- arglist
                arglist1$data <- data1
                res0 <- try( do.call( what=lsem.estimate, args=arglist1,
                                envir=environment() ), silent=TRUE )
                if ( ! inherits(res0,'try-error') ){
                    iterate <- FALSE
                    res0_out <- list(est=res0$parameters$est,
                                        M=res0$parameters_summary$M,
                                        SD=res0$parameters_summary$SD,
                                        MAD=res0$parameters_summary$MAD,
                                        lin_slo=res0$parameters_summary$lin_slo)
                }
            }
            return(res0_out)
        }


        res_all <- sirt_parlapply(cl=cl, X=1:B, FUN=fun_lsem_permutation_test,
                        verbose=verbose, arglist=arglist)

        parallel::stopCluster(cl)

        for (bb in 1L:B){
            parameters_permutation[, bb] <- res_all[[bb]]$est
            parameters_summary_M[,bb] <- res_all[[bb]]$M
            parameters_summary_SD[,bb] <- res_all[[bb]]$SD
            parameters_summary_MAD[,bb] <- res_all[[bb]]$MAD
            parameters_summary_lin_slo[,bb] <- res_all[[bb]]$lin_slo
        }

        nn <- NA
    }
    #-----------------

    #*** create global test statistics
    teststat <- data.frame( par=parameters_summary$par )
    teststat$M <- parameters_summary$M
    teststat$SD <- parameters_summary$SD
    teststat$SD_p <- rowMeans( parameters_summary_SD >=parameters_summary$SD )
    teststat$MAD <- parameters_summary$MAD
    teststat$MAD_p <- rowMeans( parameters_summary_MAD >=parameters_summary$MAD )
    teststat$lin_slo <- parameters_summary$lin_slo
    p1 <- rowMeans( parameters_summary_lin_slo >=parameters_summary$lin_slo )
    p2 <- rowMeans( parameters_summary_lin_slo <=parameters_summary$lin_slo )
    teststat$lin_slo_p <- 2*ifelse( p1 < p2, p1, p2 )
    teststat$lin_slo_p <- ifelse( teststat$lin_slo_p > 1, 1, teststat$lin_slo_p )
    #********************
    # pointwise statistics
    rownames(parameters_permutation) <- rownames(parameters)

    parameters_pointwise_test <- parameters[, c('grid_index',
                                    'moderator','par','parindex' ) ]
    parindex_ <- parameters$parindex
    parameters_pointwise_test$est <- parameters$est - parameters_summary$M[ parindex_ ]
    par_pointwise_perm <- parameters_permutation -     parameters_summary_M[ parindex_, ]

    p1 <- rowMeans( parameters_pointwise_test$est >=par_pointwise_perm )
    p2 <- rowMeans( parameters_pointwise_test$est <=par_pointwise_perm )
    parameters_pointwise_test$p <- 2*ifelse( p1 < p2, p1, p2 )
    parameters_pointwise_test$p <- ifelse( parameters_pointwise_test$p > 1,
                                        1, parameters_pointwise_test$p )

    #- non-convergence rates
    nonconverged_rate <- 100*nn / (B+nn)

    s2 <- Sys.time()

    res <- list( teststat=teststat,
                    parameters_pointwise_test=parameters_pointwise_test,
                    parameters=parameters,
                    parameters_permutation=parameters_permutation,
                    parameters_summary=parameters_summary,
                    parameters_summary_M=parameters_summary_M,
                    parameters_summary_SD=parameters_summary_SD,
                    parameters_summary_MAD=parameters_summary_MAD,
                    parameters_summary_lin_slo=parameters_summary_lin_slo,
                    par_pointwise_perm=par_pointwise_perm,
                    moderator.density=object$moderator.density,
                    moderator=object$moderator,
                    moderator.grid=object$moderator.grid,
                    h=object$h, bw=object$bw, N=object$N,
                    nonconverged_rate=nonconverged_rate,
                    use_lavaan_survey=use_lavaan_survey,
                    B=B, s1=s1, s2=s2, lavmodel=object$lavmodel,
                    n.core=n.core, cl.type=cl.type, CALL=CALL )
    class(res) <- 'lsem.permutationTest'
    return(res)
}
