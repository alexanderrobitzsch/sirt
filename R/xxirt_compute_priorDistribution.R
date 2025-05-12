## File Name: xxirt_compute_priorDistribution.R
## File Version: 0.148


xxirt_compute_priorDistribution <- function( Theta, customTheta, G )
{
    P_Theta <- customTheta$P
    arg_Theta <- list( Theta=Theta, par=customTheta$par, G=G )
    if (customTheta$person_covariates){
        arg_Theta$X <- customTheta$X
    }
    prior_Theta <- do.call( what=P_Theta, args=arg_Theta )
    return(prior_Theta)
}
