## File Name: rasch_mirtlc_mstep_mlc1.R
## File Version: 0.16


#-- M-step rasch mirtlc
rasch_mirtlc_mstep_mlc1 <- function( pjk, n.k, r.jk, n.jk, G, Nclasses,
            theta.k, b, a, I, ref.item, mstep.maxit,
            des.theta, des.b, theta.fixed, theta.normal, f.qk.yi,D,
            distribution.trait, est.a, Qmatrix,modeltype, range.b, range.a,
            iter, fac.iter)
{
    if (G==1){
        pi.k <- n.k / sum(n.k)
    }
    if ( G> 1){
        pi.k <- n.k / matrix( colSums(n.k), nrow=Nclasses, ncol=G, byrow=TRUE)
    }

    #*** perform logistic regression
    if (G==1){
        c1 <- r.jk[,,1]
        i1 <- n.jk[,,1] - c1
    }
    if (G > 1){
        c1 <- apply( r.jk, c(1,2), sum )
        i1 <- apply( n.jk - r.jk, c(1,2), sum )
    }

    # class weights correct response
    wc1 <- matrix( c1, nrow=I*Nclasses, 1 )
    wi1 <- matrix( i1, nrow=I*Nclasses, 1 )
    # outcome
    y <- c( rep(1,I*Nclasses), rep(0,I*Nclasses) )
    wgt <- c( wc1, wi1 )

    if ( theta.fixed ){
        if (D==1){
            theta.offset <- rep( theta.k, each=I )
            theta.offset <- c( theta.offset, theta.offset )
        }
        if (D>1){
            tk1 <- matrix( theta.k, nrow=1, ncol=ncol(des.theta) )
            tk1 <- matrix( tk1, nrow=nrow(des.theta), ncol=ncol(tk1), byrow=TRUE)
            theta.offset <- des.theta * tk1
            theta.offset <- rowSums( theta.offset )
        }
    }
    if ( ( ! theta.fixed ) ){
        b_old <- b
        theta_old <- theta.k
        # structure theta
        # theta_1: class_1, ..., class_D
        # theta_2: class_1, ..., class_D
        # ... theta_D: class_1, ..., class_D
        # unconstrained theta optimization
        if (modeltype=="MLC2"){
            des1 <- matrix( a, nrow=nrow(des.theta), ncol=1)
            des.theta <- des.theta * outer( des1[,1], rep(1,ncol(des.theta) ) )
            des.b <- des.b * outer( des1[,1], rep(1,ncol(des.b) ) )
        }
            mod2 <- stats::glm( y ~ 0 + des.theta + des.b, weights=wgt, family="binomial",
                        control=list(maxit=mstep.maxit ) )
        if (D==1){
            theta.k <- coef(mod2)[ 1:Nclasses ]    # theta
            b0 <- coef(mod2)[ -c( 1:Nclasses ) ]
        }
        if (D>1){
            theta.k0 <- theta.k
            theta.k <- coef(mod2)[ 1:(D*Nclasses) ]    # theta
            theta.k <- matrix( theta.k, nrow=Nclasses, ncol=D )
            theta.change <- theta.k - theta_old
            increment <- .9^( iter^fac.iter)
            theta.change <- ifelse( abs( theta.change ) > increment, sign(theta.change)*increment, theta.change )
            theta.k <- theta_old + theta.change
            b0 <- coef(mod2)[ -c( 1:(D*Nclasses) ) ]
        }
        b[ setdiff( 1:I, ref.item ) ] <- b0
        b.change <- b - b_old
        b.increment <- .9^( iter^fac.iter)
        b.change <- ifelse( abs( b.change ) > b.increment, sign(b.change)*b.increment, b.change )
        b <- b_old + b.change
    }
    if ( theta.fixed ){
        # constrained theta optimization
        if (modeltype=="MLC2"){
            des1 <- matrix( a, nrow=nrow(des.theta), ncol=1)
            theta.offset <- des1[,1] * theta.offset
            des.b <- des.b * outer( des1[,1], rep(1,ncol(des.b) ) )
        }
        mod2 <- stats::glm( y ~ 0 + offset(theta.offset) + des.b, weights=wgt, family="binomial",
                        control=list(maxit=mstep.maxit ) )
        b0 <- coef(mod2)
        b[ setdiff( 1:I, ref.item ) ] <- b0

        #***** normal distribution assumption
        if (distribution.trait=="normal" & ( D==1) ){ # D=1
            for (gg in 1:G){
                pik1 <-    pi.k[,gg]
                m1 <- sum( theta.k * pik1 )
                sd1 <- sqrt( sum( theta.k^2 * pik1 ) - m1^2 )
                pi.k[,gg] <- sirt_dnorm_discrete( theta.k, mean=m1, sd=sd1 )
            }
        }

        #--- log linear smoothing
        if (distribution.trait %in% c("smooth2", "smooth3", "smooth4") & ( D==1) ){ # D=1
            for (gg in 1:G){
                pik1 <-    pi.k[,gg]
                pik1 <- pik1 + 1e-10
                lpik1 <- log( pik1 )
                tk <- theta.k
                if ( distribution.trait=="smooth2"){
                    formula1 <- lpik1 ~ tk + I(tk^2)
                }
                if ( distribution.trait=="smooth3"){
                    formula1 <- lpik1 ~ tk + I(tk^2) + I(tk^3)
                }
                if ( distribution.trait=="smooth4"){
                    formula1 <- lpik1 ~ tk + I(tk^2) + I(tk^3) + I(tk^4)
                }
                mod <- stats::lm( formula1, weights=pik1 )
                pik2 <- exp( fitted(mod))
                pi.k[,gg] <- pik2 / sum(pik2)
            }
        }
    }

    # range restrictions
    b[ b < range.b[1] ] <- range.b[1]
    b[ b > range.b[2] ] <- range.b[2]
    ## 2PL estimation
    if (modeltype=="MLC2" ){
        a <- rasch_mirtlc_est_a( theta.k=theta.k, b=b, fixed.a=a,
                pjk=pjk, alpha1=0, alpha2=0, h=1e-4, G=G, I=I, r.jk=r.jk,
                n.jk=n.jk, est.a=est.a, Qmatrix=Qmatrix,
                iter=iter, fac.iter=fac.iter, dimensions=dimensions )
        a[ a < range.a[1] ] <- range.a[1]
        a[ a > range.a[2] ] <- range.a[2]
    }

    #-- output
    res <- list( pi.k=pi.k, pjk=pjk, theta.k=theta.k, b=b, a=a )
    return(res)
}


.m.step.mirtlc.mlc1 <- rasch_mirtlc_mstep_mlc1
