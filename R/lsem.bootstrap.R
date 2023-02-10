## File Name: lsem.bootstrap.R
## File Version: 0.333


lsem.bootstrap <- function(object, R=100, verbose=TRUE, cluster=NULL,
        repl_design=NULL, repl_factor=NULL, use_starting_values=TRUE)
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
    lsem_args$se <- "none"
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

    #-- loop over bootstrap samples
    lsem_bootstrap_print_start(verbose=verbose)

    repl_design_used <- matrix(NA, nrow=nrow(data), ncol=R)

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
        if ( ! inherits(mod1,"try-error" ) ){
            parameters_boot[,rr] <- mod1$parameters$est
            fitstats_joint_boot[,rr] <- mod1$fitstats_joint$value
            parameters_var_boot[,rr] <- mod1$parameters_summary$SD^2
            repl_design_used[,rr] <- res$repl_vector
            rr <- rr + 1
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
    object$R <- R
    object$class_boot <- TRUE
    object$fitstats_joint <- fitstats_joint
    object$repl_design <- repl_design
    object$repl_design_used <- repl_design_used
    s2 <- Sys.time()
    object$s1 <- s1
    object$s2 <- s2
    object$time <- s2-s1

    #--- output
    class(object) <- c("lsem","lsem.boot")
    return(object)
}
