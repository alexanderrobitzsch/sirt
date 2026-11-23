## File Name: sirt_dinvgamma.R
## File Version: 0.01


#--- copied from MCMCpack package
sirt_dinvgamma <- function (x, shape, scale=1)
{
    if (shape <=0 | scale <=0) {
        stop('Shape or scale parameter negative in dinvgamma().\n')
    }
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
    return(exp(log.density))
}
