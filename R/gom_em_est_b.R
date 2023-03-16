## File Name: gom_em_est_b.R
## File Version: 0.06



#--- estimation of b parameters in Rasch GOM model
gom_em_est_b <- function( lambda, I, K, n.ik, b, theta0.k, numdiff.parm=.001,
        max.increment, theta.k, msteps, mstepconv, eps=.001, progress=progress )
{
    h <- numdiff.parm
    diffindex <- 1:I
    if (progress){
        cat("  M steps b parameter |")
    }
    an.ik <- aperm( n.ik, c(2,3,1) )
    it <- 0
    conv1 <- 1000
    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        b0 <- b
        pjk <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=b,
                                    theta0.k=theta0.k )$probsL
        pjk1 <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=b+h,
                                    theta0.k=theta0.k )$probsL
        pjk2 <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=b-h,
                                    theta0.k=theta0.k)$probsL
        # numerical differentiation
        res <- gom_em_numdiff_index( pjk, pjk1, pjk2, an.ik, diffindex,
                        max.increment=max.increment, numdiff.parm )
        b <- b + res$increment
        #b <- b - mean(b)
        conv1 <- max( abs( b - b0 ) )
        it <- it+1
        if (progress){
            cat("-")
        }
    }
    if (progress){
        cat(" ", it, "Step(s) \n")
    }
    res <- list(b=b, se.b=sqrt(-1/res$d2), ll=sum(res$ll0) )
    return(res)
}


