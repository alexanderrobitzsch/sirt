## File Name: lsem_bootstrap_inference.R
## File Version: 0.08

lsem_bootstrap_inference <- function(parameters_boot, est, repl_factor=NULL)
{
    R <- ncol(parameters_boot)
    est_boot <- rowMeans(parameters_boot, na.rm=TRUE)
    if (is.null(repl_factor)){
        repl_factor <- 1/(R-1)
    }
    se_boot <- sqrt( rowSums( ( parameters_boot - est_boot )^2 ) * repl_factor )
    bias_boot <- (est_boot - est)*repl_factor*(R-1)
    est_bc <- est - bias_boot
    #-- output
    res <- list(mean_boot=est_boot, se_boot=se_boot, est_bc=est_bc,
                bias_boot=bias_boot, est=est)
    return(res)
}
