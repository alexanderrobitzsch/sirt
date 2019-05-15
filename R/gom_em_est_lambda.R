## File Name: gom_em_est_lambda.R
## File Version: 0.196


# estimation of lambda parameters
gom_em_est_lambda <- function( lambda, I, K, n.ik, numdiff.parm=.001,
        max.increment=1, theta.k, msteps, mstepconv, eps=2*numdiff.parm, progress=TRUE)
{
    h <- numdiff.parm
    diffindex <- 1:I
    Q0 <- 0*lambda
    se.lambda <- Q0
    an.ik <- aperm( n.ik, c(2,3,1) )
    if (progress){
        cat("  M steps lambda parameter |")
    }
    it <- 0
    conv1 <- 1000

    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        lambda0 <- lambda
        # nu <- stats::qlogis(lambda)
        # nu0 <- nu
        for (kk in 1:K){
            Q1 <- Q0
            Q1[,kk] <- 1
            # parm_temp <- lambda0
            parm_temp <- lambda
            # lambda00 <- stats::plogis(nu0)
            # lambda01 <- stats::plogis(nu0+h*Q1)
            # lambda02 <- stats::plogis(nu0-h*Q1)
            lambda00 <- parm_temp
            lambda01 <- parm_temp + h*Q1
            lambda02 <- parm_temp - h*Q1

            pjk <- gom_em_calc_probs( lambda=lambda00, theta.k=theta.k )$probsL
            pjk1 <- gom_em_calc_probs( lambda=lambda01, theta.k=theta.k )$probsL
            pjk2 <- gom_em_calc_probs( lambda=lambda02, theta.k=theta.k )$probsL
            res <- gom_em_numdiff_index( pjk, pjk1, pjk2, an.ik, diffindex,
                        max.increment=max.increment, numdiff.parm )
            increment <- Q1*matrix( res$increment, nrow=I, ncol=K)
            lambda <- lambda + increment
            # nu <- nu + increment
            # lambda <- stats::plogis(nu)
            lambda[ lambda[,kk] < eps, kk ] <- eps
            lambda[ lambda[,kk] > 1-eps, kk ] <- 1 - eps
            se.lambda[,kk] <- sqrt(-1/res$d2)
        }
        conv1 <- max( abs( lambda - lambda0 ) )
        it <- it+1
        if (progress){
            cat("-")
        }
    }
    if (progress){
        cat(" ", it, "Step(s) \n")
    }
    res <- list(lambda=lambda, se.lambda=se.lambda, ll=sum(res$ll0) )
    return(res)
}
