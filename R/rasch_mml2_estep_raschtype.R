## File Name: rasch_mml2_estep_raschtype.R
## File Version: 0.06


#--- E Step Raschtype Model
rasch_mml2_estep_raschtype <- function( dat1, dat2, dat2.resp,
        theta.k, pi.k, I, n, b, fixed.a, fixed.c,  fixed.d, alpha1, alpha2,
        group, pseudoll, f.qk.yi=NULL )
{
    # arrange groups
    if ( is.null(group) ){
        group <- rep( 1, nrow(dat1))
    }
    G <- length( unique(group) )
    TP <- length(theta.k)

    # probabilities of correct item at theta_k
    pjk <- .prob.raschtype.genlogis( theta.k, b, alpha1, alpha2, fixed.a )
# cat("   probs") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1

    fixed.c.M <- sirt_matrix2( x=fixed.c, nrow=nrow(pjk) )
    fixed.d.M <- sirt_matrix2( x=fixed.d, nrow=nrow(pjk) )
    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
# cat("   probs compute") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1

    if ( is.null( f.qk.yi ) ){
        pjkt <- t(pjk)
        pjkL <- array( NA, dim=c( I, 2, TP ) )
        pjkL[,1,] <- 1 - pjkt
        pjkL[,2,] <- pjkt
        probsM <- matrix( aperm( pjkL, c(2,1,3) ), nrow=I*2, ncol=TP )
        f.yi.qk <- mml_calc_like( dat2=dat2, dat2resp=dat2.resp, probs=probsM,
                            pseudoll=pseudoll )$fyiqk
#cat("   calc likelihood") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1
        # f.qk.yi <- 0 * f.yi.qk
        f.qk.yi <- matrix( 0, nrow=nrow(f.yi.qk), ncol=ncol(f.yi.qk) )
        if ( G==1 ){
            pi.k <- matrix( pi.k, ncol=1 )
        }
        for ( gg in 1:G){
            ind_gg <- which(group==gg)
            f.qk.yi[ind_gg,] <- f.yi.qk[ind_gg,]*sirt_matrix2(x=pi.k[,gg], nrow=length(ind_gg))
        }
        f.qk.yi <- f.qk.yi / rowSums(f.qk.yi)
    }

#cat("   calc posterior") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1

    # expected counts at theta.k
    n.k <- matrix( 0, length(theta.k), G )
    r.jk <- n.jk <- array( 0, dim=c( ncol(dat2), length(theta.k), G) )
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

#cat("   calc expected counts") ; vv1 <- Sys.time(); print(vv1-vv0) ; vv0 <- vv1

    #-- output
    res <- list( n.k=n.k, n.jk=n.jk, r.jk=r.jk, f.qk.yi=f.qk.yi, pjk=pjk,
                f.yi.qk=f.yi.qk, ll=sum(ll) )
    return(res)
}

.e.step.raschtype <- rasch_mml2_estep_raschtype
