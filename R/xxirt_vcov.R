
vcov.xxirt <- function( object , ...)
{
	res <- xxirt_hessian( object )
	return( solve(-res) )
}
