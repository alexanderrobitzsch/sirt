## File Name: xxirt_compute_priorDistribution.R
## File Version: 0.06

###############################################################################
xxirt_compute_priorDistribution <- function( Theta , customTheta , G )
{
	P_Theta <- customTheta$P	
	arg_Theta <- list( "Theta" = Theta , "par" = customTheta$par , "G" = G )
    prior_Theta <- do.call( P_Theta , arg_Theta )
	return(prior_Theta)
}
###############################################################################				
