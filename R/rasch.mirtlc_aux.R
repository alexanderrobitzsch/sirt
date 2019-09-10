## File Name: rasch.mirtlc_aux.R
## File Version: 91.26

#*************************************************************************************
# E Step Multidimensional Latent Class Rasch Model                                 #
.e.step.mirtlc.mlc1 <- function( dat1, dat2, dat2.resp, pi.k, pjk, I,
                   b, a, group, G, theta.k,  D, dimensions, Qmatrix,
                   f.qk.yi=NULL  ){
    #...................................
        #***
        # array notation of probabilities
        if ( D==1){
            pjk <- .prob.raschtype.genlogis( theta.k, b,
                    alpha1=0, alpha2=0, fixed.a=a )
                    }
        if (D>1){
            thetaPred <- theta.k %*% t(Qmatrix )
            bPred <- matrix( b, nrow=nrow(theta.k), ncol=I, byrow=TRUE)
            aPred <- matrix( a, nrow=nrow(theta.k), ncol=I, byrow=TRUE)
            pjk <- stats::plogis( aPred*(thetaPred - bPred ) )
                }
#        pjkL <- array( NA, dim=c(2, nrow(pjk), ncol(pjk) ) )
#        pjkL[1,,] <- 1 - pjk
#        pjkL[2,,] <- pjk
        if (D==1){ NT <- length(theta.k)  } else {NT <- nrow(theta.k) }
        TP <- NT
        pjkt <- t(pjk)
        pjkL <- array( NA, dim=c( I, 2, TP  ) )
        pjkL[,1,] <- 1 - pjkt
        pjkL[,2,] <- pjkt
        probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
        f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp,
                    probs=probsM )$fyiqk

#        f.yi.qk <- matrix( 1, nrow(dat2), NT )
#        for (ii in 1:ncol(dat2)){
        #    ii <- 1
#            ind.ii <- which( dat2.resp[,ii]==1 )
#            f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * pjkL[ dat2[ind.ii,ii]+1,,ii]
#                        }
        #******
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k, ncol=1 ) }
    if (G>1){
        for ( gg in 1:G){
            f.qk.yi[ group==gg, ] <- f.yi.qk[ group==gg, ] *
                    outer( rep( 1, nrow(dat2[ group==gg,]) ), pi.k[,gg] )
                        }
                    }
    if (G==1){
            f.qk.yi <- f.yi.qk * matrix( pi.k[,1], nrow=nrow(f.yi.qk), ncol=nrow(pi.k), byrow=TRUE )
            }
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
    # expected counts at theta.k
    if (D==1){ NT <- length(theta.k) } else { NT <- nrow(theta.k ) }
    n.k <- matrix( 0, NT, G )
    r.jk <- n.jk <- array( 0, dim=c( ncol(dat2), NT, G) )
    ll <- rep(0,G)

    for (gg in 1:G){
        ind.gg <- which( group==gg )
        res <- mml_raschtype_counts( dat2=dat2[ind.gg,], dat2resp=dat2.resp[ind.gg,],
                    dat1=dat1[ind.gg,2], fqkyi=f.qk.yi[ind.gg,],
                    pik=pi.k[,gg], fyiqk=f.yi.qk[ind.gg,]  )
        n.k[,gg] <- res$nk
        n.jk[,,gg] <- res$njk
        r.jk[,,gg] <- res$rjk
        ll[gg] <- res$ll
                }
    res <- list( "n.k"=n.k, "n.jk"=n.jk, "r.jk"=r.jk, "f.qk.yi"=f.qk.yi, "pjk"=pjk,
            "f.yi.qk"=f.yi.qk, "ll"=sum(ll) )
    return(res)
    }
#*************************************************************************************





#*************************************************************************************
# E Step Raschtype Model                                                        #
.e.step.mirtlc.lc <- function( dat1, dat2, dat2.resp, pi.k, pjk, I,
                   group, G, theta.k,  f.qk.yi=NULL  ){
    #...................................
        #***
        # array notation of probabilities
        TP <- nrow(pjk)
        pjkt <- t(pjk)
        pjkL <- array( NA, dim=c( I, 2, TP  ) )
        pjkL[,1,] <- 1 - pjkt
        pjkL[,2,] <- pjkt
        probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
        f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp,
                    probs=probsM )$fyiqk

        #******
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k, ncol=1 ) }
    if (G>1){
        for ( gg in 1:G){
            ind.gg <- group==gg
            f.qk.yi[ ind.gg, ] <- f.yi.qk[ ind.gg, ] * matrix( pi.k[,gg], nrow=length(ind.gg),
                            ncol=nrow(pi.k), byrow=TRUE )
                        }
                }
    if (G==1){
            f.qk.yi <- f.yi.qk * matrix( pi.k[,1], nrow=nrow(f.yi.qk), ncol=nrow(pi.k), byrow=TRUE )
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
#        n.jk[,,gg] <- ( t(dat2.resp[ind.gg,]) * outer( rep(1,I), dat1[ind.gg,2] ) ) %*%
#                    f.qk.yi[ ind.gg, ]
        f.qk.yi.gg <- f.qk.yi[ ind.gg, ]

        dat1.gg <- matrix( dat1[ind.gg,2], nrow=I, ncol=sum(ind.gg), byrow=TRUE )
        M1 <- ( t(dat2.resp[ind.gg,]) * dat1.gg )
        n.jk[,,gg] <-  M1 %*% f.qk.yi.gg

        # compute r.jk (expected counts for correct item responses at theta.k for item j
        r.jk[,,gg] <- ( t(dat2[ind.gg,]) * M1 ) %*% f.qk.yi.gg
        # compute log-Likelihood
        ll[gg] <- sum( dat1[ind.gg,2] * log( rowSums( f.yi.qk[ind.gg,] *
                    outer( rep( 1,nrow(f.yi.qk[ind.gg,,drop=FALSE]) ), pi.k[,gg] ) ) ) )
                }
    res <- list( "n.k"=n.k, "n.jk"=n.jk, "r.jk"=r.jk, "f.qk.yi"=f.qk.yi, "pjk"=pjk,
            "f.yi.qk"=f.yi.qk, "ll"=sum(ll) )
    return(res)
    }
#*************************************************************************************




########################################################
# calculate class probabilities
.m.step.mirtlc.lc <- function( pjk, n.k, r.jk, n.jk, G, Nclasses ){
    if (G==1){
        pi.k <- n.k / sum( n.k )
        pi.k <- matrix( pi.k, nrow=ncol(pjk), ncol=nrow(pjk) )
            }
    if ( G> 1){
        pi.k <- n.k / matrix( colSums(n.k ), nrow=Nclasses, ncol=G, byrow=TRUE)
            }
    for (cc in 1:Nclasses ){
    # cc <- 1
    if (G==1){
            pjk[ cc, ] <- r.jk[, cc, 1] / n.jk[, cc, 1]
                }
    if (G>1){
            pjk[ cc, ] <- rowSums( r.jk[, cc, ] ) / rowSums( n.jk[, cc, ]  )
                }
                    }
    res <- list( "pi.k"=pi.k, "pjk"=pjk )
    return(res)
            }
########################################################

