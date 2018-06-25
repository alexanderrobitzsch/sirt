## File Name: rm_numdiff_trim_increment.R
## File Version: 0.07

rm_numdiff_trim_increment <- function( increment, max.increment, eps2 )
{
    ci <- ceiling( abs(increment) / ( abs( max.increment) + eps2 ) )
    increment <- ifelse( abs(increment) > abs(max.increment),
                        increment/(2*ci), increment )
    return(increment)
}
