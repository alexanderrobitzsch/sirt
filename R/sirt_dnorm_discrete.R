## File Name: sirt_dnorm_discrete.R
## File Version: 0.04
## File Last Change: 2019-05-03

sirt_dnorm_discrete <- function(x, mean=0, sd=1, ...)
{
    fx <- stats::dnorm(x, mean=mean, sd=sd, ...)
    fx <- fx / sum(fx)
    return(fx)
}
