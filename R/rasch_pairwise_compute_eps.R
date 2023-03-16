## File Name: rasch_pairwise_compute_eps.R
## File Version: 0.02
## File Last Change: 2021-03-28


rasch_pairwise_compute_eps <- function(x)
{
    I <- length(x)+1
    eps <- rep(1,I)
    eps[2:I] <- exp(-x)
    return(eps)
}
