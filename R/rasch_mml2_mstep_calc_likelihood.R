## File Name: rasch_mml2_mstep_calc_likelihood.R
## File Version: 0.02

rasch_mml2_mstep_calc_likelihood <- function(G, pjk.M, n.jk, r.jk )
{
    qjk.M <- 1 - pjk.M
    ll0 <- rep(0,G)
    for (gg in 1:G){
        ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
    }
    res <- sum(ll0)
    return(res)
}
