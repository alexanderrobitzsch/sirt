## File Name: sirt_squeeze.R
## File Version: 0.03

sirt_squeeze <- function(x, lower=NULL, upper=NULL)
{
    if (!is.null(lower)){
        x <- ifelse( x<lower, lower, x)
    }
    if (!is.null(upper)){
        x <- ifelse( x>upper, upper, x)
    }
    return(x)
}
