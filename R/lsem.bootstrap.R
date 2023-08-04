## File Name: lsem.bootstrap.R
## File Version: 0.431


lsem.bootstrap <- function(object, R=100, verbose=TRUE, cluster=NULL,
        repl_design=NULL, repl_factor=NULL, use_starting_values=TRUE,
        n.core=1, cl.type="PSOCK")
{
    # fix seed locally
    s1 <- Sys.time()

    #* do not fix seed
    if (FALSE){
        old <- .Random.seed
        on.exit({.Random.seed <<- old})
        # set.seed(seed)
    }

    if (!is.null(repl_design)){
        R <- ncol(repl_design)
    }

    lsem_args <- object$lsem_args
    lsem_args$se <- 'none'
    lsem_args$verbose <- FALSE
    lsem_args$sampling_weights <- object$sampling_weights

    #- use starting values
    if (use_starting_values){
        partable_joint <- object$partable_joint
        if ( !is.null(partable_joint)){
            partable_joint$start <- partable_joint$est
            lsem_args$partable_joint <- partable_joint
        }
    }

    sampling_weights <- object$sampling_weights
    data <- object$data
    est_joint <- object$est_joint
    fitstats_joint <- object$fitstats_joint

    #-- bootstrap parameters_summary
    parameters_summary <- object$parameters_summary
    parameters_var_boot <- matrix(NA, nrow=nrow(parameters_summary), ncol=R)

    #-- create output arguments
    parameters <- object$parameters
    NP <- nrow(parameters)
    parameters_boot <- matrix(NA, nrow=NP, ncol=R)
    rownames(parameters_boot) <- rownames(parameters)
    if (est_joint){
        fitstats_joint_boot <- matrix(NA, nrow=nrow(fitstats_joint), ncol=R)
        rownames(fitstats_joint_boot) <- paste(fitstats_joint$stat)
    } else {
        fitstats_joint_boot <- NULL
    }

    #- moderator grid
    NG <- nrow(object$moderator.density)
    moderator_density_boot <- matrix(NA, nrow=NG, ncol=R)

    #-- loop over bootstrap samples
    lsem_bootstrap_print_start(verbose=verbose)

    repl_design_used <- matrix(NA, nrow=nrow(data), ncol=R)

    if (n.core==1){
        rr <- 1
        while (rr<=R){
            lsem_bootstrap_print_progress(rr=rr, verbose=verbose, R=R)
            #- draw bootstrap sample
            res <- lsem_bootstrap_draw_bootstrap_sample(data=data,
                                sampling_weights=sampling_weights, lsem_args=lsem_args,
                                cluster=cluster, repl_design=repl_design, rr=rr)
            lsem_args1 <- res$lsem_args1
            #- fit model
            mod1 <- try( do.call(what=lsem.estimate, args=lsem_args1), silent=TRUE)

            #- output collection
            if ( ! inherits(mod1,'try-error' ) ){
                parameters_boot[,rr] <- mod1$parameters$est
                fitstats_joint_boot[,rr] <- mod1$fitstats_joint$value
                parameters_var_boot[,rr] <- mod1$parameters_summary$SD^2
                repl_design_used[,rr] <- res$repl_vector
                moderator_density_boot[,rr] <- mod1$moderator.density$wgt
                rr <- rr + 1
            }
        }
    }
    if (n.core>1){
        m <- parallel::detectCores()
        if(n.core>m-1){
            n.core <- m-1
        }
        cl <- parallel::makeCluster(spec=n.core, type=cl.type)

        arglist <- list( data=data, sampling_weights=sampling_weights,
                        lsem_args=lsem_args, cluster=cluster, repl_design=repl_design )

        varlist <- c('arglist', 'lsem_bootstrap_draw_bootstrap_sample')
        parallel::clusterExport(cl=cl, varlist=varlist, envir=environment() )
        parallel::clusterEvalQ(cl, expr=library(sirt))

        fun_lsem_bootstrap <- function(i, arglist){
            iterate <- TRUE
            data <- arglist$data
            sampling_weights <- arglist$sampling_weights
            lsem_args <- arglist$lsem_args
            cluster <- arglist$cluster
            repl_design <- arglist$repl_design
            #-- estimate LSEM
            while (iterate){
                res <- lsem_bootstrap_draw_bootstrap_sample(data=data,
                                sampling_weights=sampling_weights,
                                lsem_args=lsem_args,
                                cluster=cluster, repl_design=repl_design, rr=i)
                lsem_args1 <- res$lsem_args1

                mod1 <- try( do.call(what=lsem.estimate, args=lsem_args1),
                                silent=TRUE)
                res_out <- list()
                if ( ! inherits(mod1,'try-error' ) ){
                    res_out$parameters_boot <- mod1$parameters$est
                    res_out$fitstats_joint_boot <- mod1$fitstats_joint$value
                    res_out$parameters_var_boot <- mod1$parameters_summary$SD^2
                    res_out$moderator_density_boot <- mod1$moderator.density$wgt
                    res_out$repl_design_used <- res$repl_vector
                    iterate <- FALSE
                }
            }
            return(res_out)
        }

        res_all <- sirt_parlapply(cl=cl, X=1:R, FUN=fun_lsem_bootstrap,
                        verbose=verbose, arglist=arglist)
        parallel::stopCluster(cl)

        for (rr in 1:R){
            res_out_rr <- res_all[[rr]]
            parameters_boot[,rr] <- res_out_rr$parameters_boot
            fitstats_joint_boot[,rr] <- res_out_rr$fitstats_joint_boot
            parameters_var_boot[,rr] <- res_out_rr$parameters_var_boot
            repl_design_used[,rr] <- res_out_rr$repl_design_used
            moderator_density_boot[,rr] <- res_out_rr$moderator_density_boot
        }

    }

    #- modify output objects
    res <- lsem_bootstrap_postproc_output( parameters=parameters,
                parameters_boot=parameters_boot, fitstats_joint=fitstats_joint,
                fitstats_joint_boot=fitstats_joint_boot, est_joint=est_joint,
                repl_factor=repl_factor, parameters_summary=parameters_summary,
                parameters_var_boot=parameters_var_boot )
    parameters <- res$parameters
    fitstats_joint <- res$fitstats_joint
    object$parameters_summary <- res$parameters_summary

    #- include new objects in output
    object$parameters_boot <- parameters_boot
    object$fitstats_joint_boot <- fitstats_joint_boot
    object$parameters <- parameters
    object$moderator_density_boot <- moderator_density_boot
    object$R <- R
    object$class_boot <- TRUE
    object$fitstats_joint <- fitstats_joint
    object$repl_design <- repl_design
    object$repl_factor <- repl_factor
    object$repl_design_used <- repl_design_used
    object$n.core <- n.core
    object$cl.type <- cl.type
    s2 <- Sys.time()
    object$s1 <- s1
    object$s2 <- s2
    object$time <- s2-s1

    #--- output
    class(object) <- c('lsem','lsem.boot')
    return(object)
}
