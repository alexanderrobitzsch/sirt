## File Name: xxirt_compute_prior_Theta_from_x.R
## File Version: 0.03

xxirt_compute_prior_Theta_from_x <- function(x, em_args)
{

    # include parameters
    customTheta <- xxirt_parTheta_include_freeParameters(
                                customTheta=em_args$customTheta,
                                x=x[ em_args$parindex_Theta ])

    #*** compute prior distribution
    prior_Theta <- xxirt_compute_priorDistribution( Theta=em_args$Theta,
                                  customTheta=customTheta, G=em_args$G )

    return(prior_Theta)
}
