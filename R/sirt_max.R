## File Name: sirt_max.R
## File Version: 0.02

sirt_max <- function(x, value=0)
{
    y <- x[ ! is.na(x) ]
    if (length(y)>0){
        z <- max(x, na.rm=TRUE)
    } else {
        z <- 0
    }
    return(z)
}
