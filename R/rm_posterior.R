## File Name: rm_posterior.R
## File Version: 0.17


#######################################################
# calculate posterior and counts
rm_posterior <- function( dat2, dat2.resp, TP, pi.k,
    K, I, probs, dat2.ind.resp )
{
    #--- calculate likelihood
    probsM <- matrix( aperm( probs, c(2,1,3) ), nrow=I*(K+1), ncol=TP )
    f.yi.qk <- rm_calclike( dat2=dat2, dat2resp=dat2.resp,    probs=probsM, K=K)$fyiqk

    #--- calculate posterior and expected counts
    prior <- matrix( pi.k, nrow=nrow(dat2), ncol=TP, byrow=TRUE )
    f.qk.yi <- f.yi.qk * prior
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
    #--- expected counts
    n.ik <- array( 0, dim=c(TP, I, K+1 ) )
    N.ik <- array( 0, dim=c(TP, I ) )
    for (kk in 1:(K+1) ){
        n.ik[,,kk] <- crossprod( f.qk.yi, dat2.ind.resp[,,kk] )
        N.ik <- N.ik + n.ik[,,kk]
    }
    pi2 <- sirt_matrix2( x=pi.k, nrow=nrow(f.yi.qk) )
    ll <- sum( log( rowSums( f.yi.qk * pi2 ) ) )

    #--- compute pi.k
    pi.k <- colMeans( f.qk.yi )

    #--- output
    res <- list( f.yi.qk=f.yi.qk, f.qk.yi=f.qk.yi,
            n.ik=n.ik, N.ik=N.ik, pi.k=pi.k, ll=ll)
    return(res)
}


.rm.posterior <- rm_posterior
