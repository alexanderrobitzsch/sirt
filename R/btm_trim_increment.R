## File Name: btm_trim_increment.R
## File Version: 0.04
## File Last Change: 2018-12-30

btm_trim_increment <- function(incr, maxincr )
{
    res <- ifelse( abs(incr) > maxincr, maxincr*sign(incr), incr )
    return(res)
}
