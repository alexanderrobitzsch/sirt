## File Name: sirt_colMaxs.R
## File Version: 0.01
## File Last Change: 2017-10-02 14:57:37

sirt_colMaxs <- function(x)
{
	return( apply( x , 2 , max , na.rm=TRUE ) )
}
