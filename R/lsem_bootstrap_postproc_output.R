## File Name: lsem_bootstrap_postproc_output.R
## File Version: 0.098
## File Last Change: 2023-03-15

lsem_bootstrap_postproc_output <- function(parameters, parameters_boot,
    fitstats_joint, fitstats_joint_boot, est_joint=FALSE, repl_factor=NULL,
    parameters_summary, parameters_var_boot)
{
    #* parameters
    res <- lsem_bootstrap_inference(parameters_boot=parameters_boot, est=parameters$est,
                repl_factor=repl_factor, bc_square=NULL)
    parameters$est_bc <- res$est_bc
    parameters$se <- res$se_boot
    parameters$z <- parameters$est / parameters$se
    parameters$pvalue <- 2*stats::pnorm(abs(-parameters$z))
    quant <- stats::qnorm(.975)
    parameters$ci.lower <- parameters$est - quant * parameters$se
    parameters$ci.upper <- parameters$est + quant * parameters$se

    #* fitstats_joint: bootstrap inference for fit statistics
    if (est_joint){
        fjb <- fitstats_joint_boot
        # bc_square <- which( rownames(fjb) %in% c('rmsea','srmr') )
        bc_square <- which( rownames(fjb) %in% c('srmr') )
        val <- fitstats_joint$value
        res <- lsem_bootstrap_inference(parameters_boot=fjb,
                        est=val, repl_factor=repl_factor,
                        bc_square=bc_square)
        fitstats_joint$value_bc <- res$est_bc
        fitstats_joint$se <- res$se_boot
    }

    #* adapt parameters_summary
    vt1 <- rowMeans(parameters_var_boot)
    vt0 <- parameters_summary$SD^2
    parameters_summary$SD_se <- apply(sqrt(parameters_var_boot), 1, stats::sd)
    w <- vt0 - (vt1-vt0)
    parameters_summary$SD_bc <- ifelse( w < 0, 0, sqrt(w) )
    h <- parameters_summary$SD_bc / (parameters_summary$SD_se+1e-100)
    h <- ifelse( abs(parameters_summary$SD) < 1e-5*abs(parameters_summary$M), 0, h )
    parameters_summary$SD_t <- h

    parameters_summary <- move_variables_df(x=parameters_summary,
                                after_var='SD', move_vars=c('SD_bc', 'SD_se', 'SD_t'))

    #-- output
    res <- list(parameters=parameters, fitstats_joint=fitstats_joint,
                    parameters_summary=parameters_summary)
    return(res)
}
