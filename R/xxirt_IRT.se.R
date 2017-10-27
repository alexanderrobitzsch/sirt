## File Name: xxirt_IRT.se.R
## File Version: 0.06


IRT.se.xxirt <- function( object , ...)
{
	c1 <- coef(object)
	v1 <- vcov(object)
	par1 <- xxirt_partable_extract_freeParameters( object$partable )
	par2 <- xxirt_parTheta_extract_freeParameters( object$customTheta )
	N1 <- length(par1)
	N2 <- length(par2)
	dfr <- data.frame("partype"= c( rep("item",N1), rep("Theta",N2) ) )
	dfr$parlabel <- names(c1)
	dfr$value <- c1
	dfr$se <- sqrt( diag(v1) )
	return(dfr)
}
