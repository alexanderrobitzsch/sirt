## File Name: dexppow.R
## File Version: 0.052


## copied from normalp::dnormp
dexppow <- function (x, mu=0, sigmap=1, pow=2, log=FALSE)
{
    p <- pow
    cost <- 2 * p^(1/p) * gamma(1 + 1/p) * sigmap
    expon1 <- (abs(x - mu))^p
    expon2 <- p * sigmap^p
    dsty <- (1/cost) * exp(-expon1/expon2)
    if (log){
        dsty <- log(dsty)
    }
    return(dsty)
}
