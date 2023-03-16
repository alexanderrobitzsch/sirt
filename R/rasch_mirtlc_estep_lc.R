## File Name: rasch_mirtlc_estep_lc.R
## File Version: 0.24
## File Last Change: 2019-09-14


#--- E Step rasch_mirtlc
rasch_mirtlc_estep_lc <- function( dat1, dat2, dat2.resp, pi.k, pjk, I,
        group, G, theta.k,  f.qk.yi=NULL  )
{
    #--- array notation of probabilities
    TP <- nrow(pjk)
    pjkt <- t(pjk)
    pjkL <- array( NA, dim=c( I, 2, TP ))
    pjkL[,1,] <- 1 - pjkt
    pjkL[,2,] <- pjkt
    probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
    f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp, probs=probsM )$fyiqk
    f.qk.yi <- 0 * f.yi.qk
    if (G==1){
        pi.k <- matrix( pi.k, ncol=1 )
    }
    if (G>1){
        for (gg in 1:G){
            ind.gg <- group==gg
            ngg <- length(ind.gg)
            pikM <- sirt_matrix2(x=pi.k[,gg], nrow=ngg)
            f.qk.yi[ ind.gg, ] <- f.yi.qk[ ind.gg, ] * pikM
        }
    }
    if (G==1){
        pikM <- sirt_matrix2( x=pi.k[,1], nrow=nrow(f.yi.qk))
        f.qk.yi <- f.yi.qk * pikM
    }
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
    # expected counts at theta.k
    n.k <- matrix( 0, length(theta.k), G )
    r.jk <- n.jk <- array( 0, dim=c( ncol(dat2), length(theta.k), G) )
    ll <- rep(0,G)
    for (gg in 1:G){
        ind.gg <- group==gg
        n.k[,gg] <- colSums( dat1[ind.gg,2] * f.qk.yi[ind.gg,,drop=FALSE]  )
        # expected counts at theta.k and item j
        f.qk.yi.gg <- f.qk.yi[ ind.gg, ]
        dat1.gg <- matrix( dat1[ind.gg,2], nrow=I, ncol=sum(ind.gg), byrow=TRUE )
        M1 <- ( t(dat2.resp[ind.gg,]) * dat1.gg )
        n.jk[,,gg] <-  M1 %*% f.qk.yi.gg
        # compute r.jk (expected counts for correct item responses at theta.k for item j
        r.jk[,,gg] <- ( t(dat2[ind.gg,]) * M1 ) %*% f.qk.yi.gg
        # compute log-Likelihood
        ngg <- length(ind.gg)
        pikM <- sirt_matrix2(x=pi.k[,gg], nrow=ngg)
        ll[gg] <- sum( dat1[ind.gg,2] * log( rowSums( f.yi.qk[ind.gg,]*pikM)))
    }
    res <- list( n.k=n.k, n.jk=n.jk, r.jk=r.jk, f.qk.yi=f.qk.yi, pjk=pjk,
                f.yi.qk=f.yi.qk, ll=sum(ll) )
    return(res)
}


.e.step.mirtlc.lc <- rasch_mirtlc_estep_lc
