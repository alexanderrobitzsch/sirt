## File Name: rm_eap_reliability.R
## File Version: 0.01

rm_eap_reliability <- function( EAP , SE_EAP )
{
	EAP.rel <- 1 - mean( SE_EAP^2 ) / ( mean( SE_EAP^2 ) + stats::var( EAP ) )	
	return(EAP.rel)
}
