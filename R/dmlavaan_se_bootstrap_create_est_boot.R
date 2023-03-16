## File Name: dmlavaan_se_bootstrap_create_est_boot.R
## File Version: 0.05


dmlavaan_se_bootstrap_create_est_boot <- function(est_boot1, est_boot2)
{
    #-- est_boot
    est_boot1 <- dmlavaan_add_suffix_column_names(x=est_boot1, suffix='_mod1')
    est_boot2 <- dmlavaan_add_suffix_column_names(x=est_boot2, suffix='_mod2')
    est_boot <- data.frame( est_boot1, est_boot2)
    colnames(est_boot) <- c(colnames(est_boot1), colnames(est_boot2) )

    #-- compute covariance matrix
    V <- stats::cov(est_boot)

    #--- output
    res <- list(est_boot=est_boot, V=V)
    return(res)
}
