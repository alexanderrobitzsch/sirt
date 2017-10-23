## File Name: rm_eap_reliability.R
## File Version: 0.01
## File Last Change: 2017-10-02 16:09:19

rm_eap_reliability <- function( EAP , SE_EAP )
{
	EAP.rel <- 1 - mean( SE_EAP^2 ) / ( mean( SE_EAP^2 ) + stats::var( EAP ) )	
	return(EAP.rel)
}
