## File Name: lsem_bootstrap_inference.R
## File Version: 0.04

lsem_bootstrap_inference <- function(parameters_boot, est)
{
    R <- ncol(parameters_boot)
    est_boot <- rowMeans(parameters_boot, na.rm=TRUE)
    se_boot <- sqrt( rowSums( ( parameters_boot - est_boot )^2 ) / (R-1) )
    bias_boot <- est_boot - est
    est_bc <- est - bias_boot
    #-- output
    res <- list(mean_boot=est_boot, se_boot=se_boot, est_bc=est_bc,
                bias_boot=bias_boot, est=est)
    return(res)
}
