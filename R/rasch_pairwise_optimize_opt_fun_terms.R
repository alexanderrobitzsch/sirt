## File Name: rasch_pairwise_optimize_opt_fun_terms.R
## File Version: 0.03
## File Last Change: 2021-03-28

rasch_pairwise_optimize_opt_fun_terms <- function(eps_horiz, eps_vert, n.ij, n.ji)
{
    eps <- eps_horiz
    I <- length(eps)
    epsM <- matrix(eps_vert, nrow=I, ncol=I, byrow=TRUE)
    t1 <- ( n.ij*eps - n.ji*epsM )^2
    t2 <- (n.ij+n.ji)*eps*epsM + 1e-7
    t3 <- t1/t2
    return(t3)
}
