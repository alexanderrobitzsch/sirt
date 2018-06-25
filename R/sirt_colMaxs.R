## File Name: sirt_colMaxs.R
## File Version: 0.02

sirt_colMaxs <- function(x)
{
    return( apply( x, 2, max, na.rm=TRUE ) )
}
