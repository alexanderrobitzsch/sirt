## File Name: rm_numdiff_trim_increment.R
## File Version: 0.01
## File Last Change: 2017-10-02 12:15:33

rm_numdiff_trim_increment <- function( increment, max.increment, eps2 )
{
	ci <- ceiling( abs(increment) / ( abs( max.increment) + eps2 ) )
    increment <- ifelse( abs(increment) > abs(max.increment)  , 
                                 increment/(2*ci) , increment )
	return(increment)
}
