## File Name: xxirt_parTheta_extract_freeParameters.R
## File Version: 0.05
## File Last Change: 2017-06-15 18:49:33


xxirt_parTheta_extract_freeParameters <- function( customTheta )
{		
	ind <- customTheta$est
	p1 <- customTheta$par[ ind ]
	names(p1) <- names(customTheta$par)[ind]
	return(p1)
}
