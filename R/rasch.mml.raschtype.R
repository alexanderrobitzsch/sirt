## File Name: rasch.mml.raschtype.R
## File Version: 2.63



#*******************************************************
# utility function squeeze
squeeze.mml2 <- function( v1, rgvec ){
    v1 <- ifelse( v1 < rgvec[1], rgvec[1], v1 )
    v1 <- ifelse( v1 > rgvec[2], rgvec[2], v1 )
    return(v1)
        }
#*******************************************************



#*************************************************************************************
# E Step Raschtype Model: multidimensional version
.e.step.raschtype.mirt <- function( dat1, dat2, dat2.resp, theta.k, pi.k, I,
                n, b, fixed.a, fixed.c,  fixed.d,
                alpha1, alpha2, group,  mu,  Sigma.cov, Qmatrix, pseudoll ){
    #...................................
    # arrange groups
# aa0 <- Sys.time()
    if ( is.null(group) ){ group <- rep( 1, nrow(dat1)) }
    G <- length( unique( group) )
    # probabilities of correct item at theta_k
    pjk <- .prob.raschtype.genlogis( theta.k, b, alpha1, alpha2, fixed.a, Qmatrix )
    fixed.c.M <- outer( rep(1,nrow(pjk)), fixed.c )
    fixed.d.M <- outer( rep(1,nrow(pjk)), fixed.d )
    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
    TP <- dim(pjk)[1]
    #***
    # array notation of probabilities
        pjkt <- t(pjk)
        pjkL <- array( NA, dim=c( I, 2, TP  ) )
        pjkL[,1,] <- 1 - pjkt
        pjkL[,2,] <- pjkt
        probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
# cat("- probs") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1
        f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp,
                    probs=probsM, pseudoll=pseudoll)$fyiqk
# cat("- likelihood") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1

    #******
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){
            pi.k <- matrix( pi.k, ncol=1 )
                    }
    for ( gg in 1:G){
        f.qk.yi[ group==gg, ] <- f.yi.qk[ group==gg, ] * outer( rep( 1, nrow(dat2[ group==gg,]) ), pi.k[,gg] )
                    }
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )
# cat("- posterior") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1
    # expected counts at theta.k
    n.k <- matrix( 0, nrow(theta.k), G )
    r.jk <- n.jk <- array( 0, dim=c( ncol(dat2), nrow(theta.k), G) )
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
# cat("- counts") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1
    res <- list( "n.k"=n.k, "n.jk"=n.jk, "r.jk"=r.jk, "f.qk.yi"=f.qk.yi, "pjk"=pjk,
            "f.yi.qk"=f.yi.qk, "ll"=sum(ll) )
    return(res)
    }
#*************************************************************************************



##################################################
# estimation of group means in the 1dim IRT model
.est.mean <- function( dat1.gg, f.yi.qk.gg, X1, pi.k, pi.k0, gg,
                            mean.trait, sd.trait, theta.k, h)
{
    ll0 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k[,gg] ) ) ) )
    pi.k2 <- pi.k0
    pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg]+h, sd=sd.trait[gg]  )
    ll1 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ))
    pi.k2 <- pi.k0
    pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg] -h, sd=sd.trait[gg]  )
    ll2 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    d2[ abs(d2) < 10^(-15) ]  <- 10^(-15)
    d.change <- - d1 / d2
    d.change <- ifelse( abs(d.change) > .1, .1*sign(d.change), d.change )
    return(d.change)
}
##################################################



##################################################
# estimation of group SD's in the 1dim IRT model
.est.sd <- function( dat1.gg, f.yi.qk.gg, X1, pi.k, pi.k0, gg,
                            mean.trait, sd.trait, theta.k, h){
        ll0 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k[,gg] ) ) ) )
        pi.k2 <- pi.k
        pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg], sd=sd.trait[gg] + h )
        ll1 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
        pi.k2 <- pi.k
        pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg], sd=sd.trait[gg] - h )
        pi.k2[,gg] <- pi.k2[,gg] / sum( pi.k2[,gg] )
        ll2 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
        d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
        d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
        d2[ abs(d2) < 10^(-15) ]  <- 10^(-15)
        d.change <- - d1 / d2
        d.change <- ifelse( abs(d.change) > .05, .05*sign(d.change), d.change )
        return(d.change )
                            }
##################################################
