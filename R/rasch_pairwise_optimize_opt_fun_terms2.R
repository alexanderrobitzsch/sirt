## File Name: rasch_pairwise_optimize_opt_fun_terms2.R
## File Version: 0.12

rasch_pairwise_optimize_opt_fun_terms2 <- function(eps_horiz, eps_vert, y.ij, y.ji,
        estimator="MINCHI")
{
    tol <- 1e-10
    eps1 <- sqrt(eps_horiz + tol)
    I <- length(eps1)
    if (estimator=="ULS"){
        eps_v <- eps_vert
    }
    if (estimator=="MINCHI"){
        eps_v <- sqrt(eps_vert+tol)
    }
    epsM1 <- matrix(eps_v, nrow=I, ncol=I, byrow=TRUE)
    if (estimator=="MINCHI"){
        h1 <- eps1/epsM1
        t3 <- ( y.ij*h1 - y.ji/h1 )^2
    }
    if (estimator=="ULS"){
        t3 <- ( y.ij/epsM1 - y.ji/eps_horiz )^2
    }
    return(t3)
}
