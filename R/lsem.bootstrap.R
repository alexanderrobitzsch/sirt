## File Name: lsem.bootstrap.R
## File Version: 0.304


lsem.bootstrap <- function(object, R=100, verbose=TRUE, cluster=NULL)
{
    lsem_args <- object$lsem_args
    lsem_args$se <- "none"
    lsem_args$verbose <- FALSE
    lsem_args$sampling_weights <- object$sampling_weights
    sampling_weights <- object$sampling_weights
    data <- object$data
    est_joint <- object$est_joint
    fitstats_joint <- object$fitstats_joint

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
    rr <- 1
    while (rr<=R){
        lsem_bootstrap_print_progress(rr=rr, verbose=verbose, R=R)
        #- draw bootstrap sample
        lsem_args1 <- lsem_bootstrap_draw_bootstrap_sample(data=data,
                            sampling_weights=sampling_weights, lsem_args=lsem_args,
                            cluster=cluster)
        #- fit model
        mod1 <- try( do.call(what=lsem.estimate, args=lsem_args1), silent=TRUE)
        #- output collection
        if ( class(mod1) !="try-error" ){
            parameters_boot[,rr] <- mod1$parameters$est
            fitstats_joint_boot[,rr] <- mod1$fitstats_joint$value
            rr <- rr + 1
        }
    }

    #- modify output objects
    res <- lsem_bootstrap_postproc_output( parameters=parameters,
                parameters_boot=parameters_boot, fitstats_joint=fitstats_joint,
                fitstats_joint_boot=fitstats_joint_boot, est_joint=est_joint )
    parameters <- res$parameters
    fitstats_joint <- res$fitstats_joint

    #- include new objects in output
    object$parameters_boot <- parameters_boot
    object$parameters <- parameters
    object$R <- R
    object$class_boot <- TRUE
    object$fitstats_joint <- fitstats_joint

    #--- output
    class(object) <- c("lsem","lsem.boot")
    return(object)
}
