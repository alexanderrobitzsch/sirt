## File Name: rasch_mirtlc_estep_mlc1.R
## File Version: 0.22
## File Last Change: 2019-10-27


# E Step Multidimensional Latent Class Rasch Model
rasch_mirtlc_estep_mlc1  <- function( dat1, dat2, dat2.resp, pi.k, pjk, I,
            b, a, group, G, theta.k,  D, dimensions, Qmatrix, f.qk.yi=NULL  )
{
    #--- array notation of probabilities
    if ( D==1){
        pjk <- prob_raschtype_genlogis( theta=theta.k, b=b,    alpha1=0,
                    alpha2=0, fixed.a=a )
    }
    if (D>1){
        thetaPred <- tcrossprod( theta.k, Qmatrix)
        TP <- nrow(theta.k)
        bPred <- sirt_matrix2(x=b, nrow=TP)
        aPred <- sirt_matrix2(x=a, nrow=TP)
        pjk <- stats::plogis( aPred*(thetaPred - bPred ) )
    }
    if (D==1){
        NT <- length(theta.k)
    } else {
        NT <- nrow(theta.k)
    }
    TP <- NT
    pjkt <- t(pjk)
    pjkL <- array( NA, dim=c( I, 2, TP  ) )
    pjkL[,1,] <- 1 - pjkt
    pjkL[,2,] <- pjkt
    probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
    f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp, probs=probsM )$fyiqk
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k, ncol=1 ) }
    if (G>1){
        for ( gg in 1:G){
            ind.gg <- which( group==gg )
            pikM <- sirt_matrix2(x=pi.k[,gg], nrow=length(ind.gg))
            f.qk.yi[ ind.gg, ] <- f.yi.qk[ ind.gg, ] * pikM
        }
    }
    if (G==1){
        f.qk.yi <- f.yi.qk * sirt_matrix2( pi.k[,1], nrow=nrow(f.yi.qk) )
    }
    f.qk.yi <- f.qk.yi / rowSums(f.qk.yi)
    # expected counts at theta.k
    if (D==1){
        NT <- length(theta.k)
    } else {
        NT <- nrow(theta.k)
    }
    n.k <- matrix( 0, nrow=NT, ncol=G )
    r.jk <- n.jk <- array( 0, dim=c( ncol(dat2), NT, G) )
    ll <- rep(0,G)

    for (gg in 1:G){
        ind.gg <- which(group==gg)
        res <- mml_raschtype_counts( dat2=dat2[ind.gg,], dat2resp=dat2.resp[ind.gg,],
                    dat1=dat1[ind.gg,2], fqkyi=f.qk.yi[ind.gg,],
                    pik=pi.k[,gg], fyiqk=f.yi.qk[ind.gg,]  )
        n.k[,gg] <- res$nk
        n.jk[,,gg] <- res$njk
        r.jk[,,gg] <- res$rjk
        ll[gg] <- res$ll
    }
    res <- list( n.k=n.k, n.jk=n.jk, r.jk=r.jk, f.qk.yi=f.qk.yi, pjk=pjk,
                    f.yi.qk=f.yi.qk, ll=sum(ll) )
    return(res)
}


.e.step.mirtlc.mlc1 <- rasch_mirtlc_estep_mlc1
