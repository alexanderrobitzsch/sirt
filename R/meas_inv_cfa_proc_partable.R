## File Name: meas_inv_cfa_proc_partable.R
## File Version: 0.05


meas_inv_cfa_proc_partable <- function(partable, items)
{
    G <- max(partable$group)
    I <- length(items)
    mu <- partable[ partable$lhs=='F' & partable$op=='~1', 'est' ]
    sigma <- sqrt(partable[ partable$lhs=='F' & partable$op=='~~', 'est' ])
    lambda <- partable[ partable$lhs=='F' & partable$op=='=~', 'est' ]
    lambda <- matrix(lambda, nrow=G, ncol=I, byrow=TRUE)
    nu <- partable[ paste(partable$rhs)=='' & partable$op=='~1' &
                        partable$lhs %in% items, 'est' ]
    nu <- matrix(nu, nrow=G, ncol=I, byrow=TRUE)
    theta <- partable[ partable$op=='~~' & partable$lhs %in% items, 'est' ]
    theta <- matrix(theta, nrow=G, ncol=I, byrow=TRUE)

    #--- output
    res <- list(mu=mu, sigma=sigma, lambda=lambda, nu=nu, theta=theta)
    return(res)
}
