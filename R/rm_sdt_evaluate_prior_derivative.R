## File Name: rm_sdt_evaluate_prior_derivative.R
## File Version: 0.06


rm_sdt_evaluate_prior_derivative <- function(partable, h)
{
    partable <- partable[ partable$est, ]
    m1 <- partable$prior_M
    sd1 <- partable$prior_SD
    y1 <- stats::dnorm(partable$value+h, mean=m1, sd=sd1, log=TRUE)
    y2 <- stats::dnorm(partable$value-h, mean=m1, sd=sd1, log=TRUE)
    res <- - ( y1 - y2 ) / (2*h)
    return(res)
}
