## File Name: xxirt_vcov.R
## File Version: 0.07

vcov.xxirt <- function( object, ...)
{
    res <- xxirt_hessian( object )
    return( solve(-res) )
}
