## File Name: rasch_mml2_estep_raschtype_mirt.R
## File Version: 1.072
## File Last Change: 2021-09-22


#*** E Step Raschtype Model: multidimensional version
rasch_mml2_estep_raschtype_mirt <- function( dat1, dat2, dat2.resp, theta.k, pi.k, I,
                n, b, fixed.a, fixed.c,  fixed.d, alpha1, alpha2, group,  mu,  Sigma.cov,
                Qmatrix, pseudoll )
{
    aa0 <- Sys.time()
    if ( is.null(group) ){
        group <- rep( 1, nrow(dat1))
    }
    G <- length( unique( group) )
    pjk <- prob_raschtype_genlogis( theta=theta.k, b=b, alpha1=alpha1,
                    alpha2=alpha2, fixed.a=fixed.a, Qmatrix=Qmatrix )
    TP <- nrow(pjk)
# cat("- probs") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1
    fixed.c.M  <- sirt_matrix2(fixed.c, nrow=TP)
    fixed.d.M  <- sirt_matrix2(fixed.d, nrow=TP)
    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk

    TP <- dim(pjk)[1]
    pjkt <- t(pjk)
    pjkL <- array( NA, dim=c( I, 2, TP))
    pjkL[,1,] <- 1 - pjkt
    pjkL[,2,] <- pjkt
    probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
    f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp,
                    probs=probsM, pseudoll=pseudoll)$fyiqk

    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){
        pi.k <- matrix( pi.k, ncol=1 )
    }
    for ( gg in 1:G){
        f.qk.yi[ group==gg, ] <- f.yi.qk[ group==gg, ] * outer( rep( 1, nrow(dat2[ group==gg,]) ), pi.k[,gg] )
    }
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )

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

    #--- output
    res <- list( n.k=n.k, n.jk=n.jk, r.jk=r.jk, f.qk.yi=f.qk.yi, pjk=pjk,
                    f.yi.qk=f.yi.qk, ll=sum(ll) )
    return(res)
}



.e.step.raschtype.mirt <- rasch_mml2_estep_raschtype_mirt


# cat("- posterior") ; aa1 <- Sys.time(); print(aa1-aa0) ; aa0 <- aa1
