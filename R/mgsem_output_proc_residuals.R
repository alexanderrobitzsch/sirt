## File Name: mgsem_output_proc_residuals.R
## File Version: 0.06


mgsem_output_proc_residuals <- function(implied, suffstat)
{
    G <- length(implied)
    residuals_groupwise <- list()
    fit_total <- list(fit_Mu=0, fit_Sigma=0)
    for (gg in 1L:G){
        N <- suffstat[[gg]]$N
        res1 <- list()
        res1$Mu <- suffstat[[gg]]$M-implied[[gg]]$Mu
        res1$Sigma <- suffstat[[gg]]$S-implied[[gg]]$Sigma

        fit_Mu <- N*sum(suffstat[[gg]]$weights_M*res1$Mu^2)
        res1$fit_Mu <- fit_Mu
        fit_total$fit_Mu <- fit_total$fit_Mu+res1$fit_Mu

        fit_Sigma <- sum(res1$Sigma^2*suffstat[[gg]]$weights_M)
        res1$fit_Sigma <- N*sum( mgsem_vech(x=fit_Sigma) )
        fit_total$fit_Sigma <- fit_total$fit_Sigma+res1$fit_Sigma

        residuals_groupwise[[gg]] <- res1
    }
    fit_total$fit_both <- fit_total$fit_Mu + fit_total$fit_Sigma

    #--- output
    res <- list(fit_total=fit_total, residuals_groupwise=residuals_groupwise)
    return(res)
}
