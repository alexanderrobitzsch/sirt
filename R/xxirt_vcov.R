## File Name: xxirt_vcov.R
## File Version: 0.09
## File Last Change: 2019-08-02

vcov.xxirt <- function( object, ...)
{
    res <- xxirt_hessian( object )
    return( solve(-res) )
}
