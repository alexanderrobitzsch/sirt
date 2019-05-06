## File Name: sirt_dnorm.R
## File Version: 0.01

sirt_dnorm <- function(x, mean=0, sd=1, ...)
{
    y <- stats::dnorm(x=x, mean=mean, sd=sd, ...)
    return(y)
}
