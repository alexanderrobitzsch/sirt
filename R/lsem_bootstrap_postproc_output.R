## File Name: lsem_bootstrap_postproc_output.R
## File Version: 0.05

lsem_bootstrap_postproc_output <- function(parameters, parameters_boot,
    fitstats_joint, fitstats_joint_boot, est_joint=FALSE, repl_factor=NULL)
{
    #* parameters
    res <- lsem_bootstrap_inference(parameters_boot=parameters_boot, est=parameters$est,
                repl_factor=repl_factor)
    parameters$est_bc <- res$est_bc
    parameters$se <- res$se_boot
    parameters$z <- parameters$est / parameters$se
    parameters$pvalue <- 2*stats::pnorm(abs(-parameters$z))
    quant <- stats::qnorm(.975)
    parameters$ci.lower <- parameters$est - quant * parameters$se
    parameters$ci.upper <- parameters$est + quant * parameters$se

    #* fitstats_joint
    if (est_joint){
        res <- lsem_bootstrap_inference(parameters_boot=fitstats_joint_boot,
                    est=fitstats_joint$value, repl_factor=repl_factor)
        fitstats_joint$value_bc <- res$est_bc
        fitstats_joint$se <- res$se_boot
    }

    #-- output
    res <- list(parameters=parameters, fitstats_joint=fitstats_joint)
    return(res)
}
