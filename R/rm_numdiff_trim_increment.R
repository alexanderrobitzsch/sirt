## File Name: rm_numdiff_trim_increment.R
## File Version: 0.12

rm_numdiff_trim_increment <- function( increment, max.increment, eps2 )
{
    aincr <- abs(increment)
    amaxincr <- abs(max.increment)
    ci <- ceiling( aincr / ( amaxincr + eps2 ) )
    increment <- ifelse( aincr > amaxincr, increment/(2*ci), increment)
    return(increment)
}
