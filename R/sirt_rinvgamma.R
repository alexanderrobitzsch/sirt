## File Name: sirt_rinvgamma.R
## File Version: 0.01

sirt_rinvgamma <- function (n, shape, scale=1)
{
    return(1/stats::rgamma(n=n, shape=shape, rate=scale))
}
