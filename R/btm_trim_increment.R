## File Name: btm_trim_increment.R
## File Version: 0.01

btm_trim_increment <- function(incr, maxincr )
{
    res <- ifelse( abs(incr) > maxincr , maxincr*sign(incr) , incr )
    return(res)
}
