## File Name: rasch_mirtlc_est_a.R
## File Version: 0.14


#-- estimation of a parameter (discrimination parameter)
rasch_mirtlc_est_a <- function( theta.k, b, fixed.a, pjk, alpha1, alpha2, h,
        G, I, r.jk, n.jk, est.a, Qmatrix, iter, fac.iter, dimensions=NULL )
{
    eps <- 000000005
    #** a
    pjk <- .prob.raschtype.genlogis( theta.k, b, alpha1, alpha2, fixed.a, Qmatrix)
    pjk <- squeeze_probs(probs=pjk, eps=eps)
    pjk.M <- t(pjk)
    qjk.M <- 1 - pjk.M
    pjk1 <- .prob.raschtype.genlogis( theta.k, b, alpha1, alpha2, fixed.a + h,Qmatrix)
    pjk1 <- squeeze_probs(probs=pjk1, eps=eps)
    pjk1.M <- t(pjk1)
    qjk1.M <- 1 - pjk1.M
    pjk2 <- .prob.raschtype.genlogis( theta.k, b, alpha1, alpha2, fixed.a - h,Qmatrix)
    pjk2 <- squeeze_probs(probs=pjk2, eps=eps)
    pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
    # first order derivative
    # f(x+h) - f(x-h)=2* f'(x) * h
    ll0 <- ll1 <- ll2 <- matrix(0, I, G)
    #-- groups
    for (gg in 1:G){
        ll0[,gg] <- rowSums( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
        ll1[,gg] <- rowSums( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
        ll2[,gg] <- rowSums( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M )     )
    }
    ll0 <- rowSums(ll0)
    ll1 <- rowSums(ll1)
    ll2 <- rowSums(ll2)
    # aggregate with respect to estimation of a
    a1 <- stats::aggregate( cbind( ll0, ll1, ll2 ), list(est.a), sum, na.rm=T)
    a1 <- a1[ a1[,1] > 0, ]
    ll0 <- a1[,2]
    ll1 <- a1[,3]
    ll2 <- a1[,4]
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    # second order derivative
    # f(x+h)+f(x-h)=2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    d2[ abs(d2) < 1e-10 ] <- 1e-10
    # change in item difficulty
    a.change <- - d1 / d2
    # dampening parameter as in tam
    old_increment <- .2^( iter^fac.iter )
    ci <- ceiling( abs(a.change) / ( abs( old_increment) + 10^(-10) ) )
    a.change <- ifelse( abs( a.change) > abs(old_increment),
                                        a.change/(2*ci), a.change )
    a.change <- a.change[ match( est.a, a1[,1] ) ]
    if ( any( est.a==0 ) ){
        a.change[ est.a==0 ] <- 0
    }
    fixed.a <- fixed.a + a.change
    fixed.a[ fixed.a < 0 ] <- 0
    return(fixed.a)
}


.mirtlc.est.a <- rasch_mirtlc_est_a
