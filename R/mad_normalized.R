## File Name: mad_normalized.R
## File Version: 0.03

mad_normalized <- function(x)
{
    x <- stats::na.omit(x)
    res <- stats::median( abs( x - stats::median(x) ) )
    res <- res/0.6745
    return(res)
}
