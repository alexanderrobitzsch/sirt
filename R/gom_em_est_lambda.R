## File Name: gom_em_est_lambda.R
## File Version: 0.206
## File Last Change: 2019-05-18


# estimation of lambda parameters
gom_em_est_lambda <- function( lambda, I, K, n.ik, numdiff.parm=.001,
        max.increment=1, theta.k, msteps, mstepconv, eps=2*numdiff.parm,
        progress=TRUE, lambda_partable)
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
    increment_mat <- d2_mat <- d1_mat <- matrix(0, nrow=I, ncol=K)

    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        lambda0 <- lambda

        for (kk in 1:K){
            Q1 <- Q0
            Q1[,kk] <- 1
            parm_temp <- lambda0
            lambda00 <- parm_temp
            lambda01 <- parm_temp + h*Q1
            lambda02 <- parm_temp - h*Q1

            pjk <- gom_em_calc_probs( lambda=lambda00, theta.k=theta.k )$probsL
            pjk1 <- gom_em_calc_probs( lambda=lambda01, theta.k=theta.k )$probsL
            pjk2 <- gom_em_calc_probs( lambda=lambda02, theta.k=theta.k )$probsL
            res <- gom_em_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, an.ik=an.ik,
                        diffindex=diffindex, max.increment=max.increment,
                        numdiff.parm=numdiff.parm )
            d1_mat[,kk] <- res$d1
            d2_mat[,kk] <- res$d2
        }
        d1a <- rowsum(x=as.vector(d1_mat), group=lambda_partable$par_index)[,1]
        d2a <- rowsum(x=as.vector(d2_mat), group=lambda_partable$par_index)[,1]
        increment <- - d1a / d2a
        increment <- sirt_trim_increment(increment=increment, max_increment=max.increment)
        increment <- matrix( increment[ lambda_partable$par_index ], nrow=I, ncol=K)
        lambda <- lambda + increment
        d2 <- matrix(d2a[ lambda_partable$par_index ], nrow=I, ncol=K)
        for (kk in 1:K){
            lambda[,kk] <- sirt_squeeze(lambda[,kk], lower=eps, upper=1-eps)
            se.lambda[,kk] <- sqrt(-1/d2[,kk])
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
