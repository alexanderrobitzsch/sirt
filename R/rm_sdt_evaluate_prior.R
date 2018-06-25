## File Name: rm_sdt_evaluate_prior.R
## File Version: 0.04

rm_sdt_evaluate_prior <- function(partable)
{
    partable <- partable[ partable$est, ]
    m1 <- partable$prior_M
    sd1 <- partable$prior_SD
    y1 <- stats::dnorm(partable$value, mean=m1, sd=sd1, log=TRUE)
    y2 <- stats::dnorm(m1, mean=m1, sd=sd1, log=TRUE)
    res <- - ( y1 - y2 )
    return(res)
}
