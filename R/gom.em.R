## File Name: gom.em.R
## File Version: 5.4012


#-- gom EM algorithm
gom.em <- function( dat, K=NULL, problevels=NULL, weights=NULL, model="GOM",
    theta0.k=seq(-5,5,len=15), xsi0.k=exp(seq(-6,3,len=15)),
    max.increment=.3, numdiff.parm=1e-4, maxdevchange=1e-6,
    globconv=1e-4, maxiter=1000, msteps=4, mstepconv=.001,
    theta_adjust=FALSE, lambda.inits=NULL, lambda.index=NULL, pi.k.inits=NULL,
    newton_raphson=TRUE, optimizer="nlminb", progress=TRUE )
{
    CALL <- match.call()
    s1 <- Sys.time()
    e1 <- environment()
    s1 <- Sys.time()
    optimization <- NULL
    if ( model=="GOMRasch"){
        K <- length(theta0.k)
        max.increment <- 1
    }
    theta_grid <- NULL
    mu <- Sigma <- b <- NULL

    #-- data preparation
    res <- gom_em_prepare_data(dat=dat, weights=weights, model=model)
    dat0 <- res$dat0
    dat2 <- res$dat2
    dat2.resp <- res$dat2.resp
    dat2.ind <- res$dat2.ind
    N <- res$N
    I <- res$I
    weights <- res$weights

    #-- index table lambda matrix
    res <- gom_em_prepare_lambda_index(lambda.index=lambda.index, I=I, K=K,
                items=colnames(dat2))
    lambda.index <- res$lambda.index
    lambda_partable <- res$lambda_partable
    lambda_index_blocks <- res$lambda_index_blocks

    # design matrices
    if (model=="GOMRasch"){
        b <- - stats::qlogis( colMeans(dat, na.rm=TRUE)  )
        theta.kM <- as.matrix( expand.grid( theta0.k, xsi0.k ))
        TP <- nrow(theta.kM)
        m1 <- exp( - ( theta.kM[,1] - matrix( theta0.k, TP, K, byrow=TRUE ) )^2 / ( 2*theta.kM[,2]^2 ) )
        theta.k <- m1 / rowSums( m1 )
        #***
        mu <- c(0, .7)
        Sigma <- as.matrix(diag(c(1, 1 )))
        pi.k <- sirt_dmvnorm_discrete( x=theta.kM, mean=mu, sigma=Sigma )
        lambda <- stats::plogis( - outer( b, theta0.k, "-" ) )
        # design matrix skill space
        #        Z <- cbind( 1, theta.kM[,1], theta.kM[,1]^2, theta.kM[,1]^3,
        #                 theta.kM[,2], theta.kM[,2]^2, theta.kM[,2]^3,
        #                 theta.kM[,1]*theta.kM[,2], theta.kM[,1]^2*theta.kM[,2])
        newton_raphson <- FALSE
    }

    if (model=="GOM"){
        theta.k <- gom_em_calc_theta(K=K, problevels=problevels)
        TP <- nrow(theta.k)
        if (is.null(pi.k.inits)){
            pi.k <- as.vector(rep(1/TP, TP ))
        } else {
            pi.k <- pi.k.inits
        }
        lambda <- gom_em_inits_lambda(I=I, K=K, lambda.inits=lambda.inits,
                        lambda_partable=lambda_partable)
        theta0.k <- NULL
    }

    if (model=="GOMnormal"){
        v1 <- list(theta0.k)
        if (K>2){
            for (kk in 2:(K-1)){
                v1[[kk]] <- theta0.k
            }
        }
        theta0 <- theta0.k
        theta_grid0 <- theta_grid <- as.matrix(as.data.frame( expand.grid(v1) ))
        TP <- nrow(theta_grid)
        mu <- rep(0,K-1)
        sigma <- rep(1,K-1)
        if (is.null(pi.k.inits)){
            pi.k <- as.vector(rep(1/TP, TP ))
        } else {
            pi.k <- pi.k.inits
        }
        theta.k <- gom_em_normal_to_membership_scores(theta_grid=theta_grid, K=K, TP=TP)
        lambda <- gom_em_inits_lambda(I=I, K=K, lambda.inits=lambda.inits,
                        lambda_partable=lambda_partable)
        theta0.k <- NULL
        fac1 <- 1
    }

    # inits
    se.b <- se.lambda <- NULL
    n.ik <- array( 0, dim=c(TP,I, 2) )
    iter <- 0
    dev0 <- dev <- 0
    conv <- devchange <- 1000
    disp <- "...........................................................\n"

    #---------- start EM algorithm --------------------
    while( ( ( maxdevchange < devchange ) & (globconv < conv) ) & ( iter < maxiter ) ){
        if (progress){
            cat(disp)
            cat("Iteration", iter+1, "   ", paste( Sys.time() ), "\n" )
        }
        # previous values
        dev0 <- dev
        pi.k0 <- pi.k
        lambda0 <- lambda
        b0 <- b
# a0 <- Sys.time()
        # calculate probabilities
        res <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=b,
                        theta0.k=theta0.k )
        probs <- res$probs
        probs <- problong2probarray( probres=probs, I=I, TP=TP )

        # calculate counts
        probsM <- matrix( aperm( probs, c(2,1,3) ), nrow=I*2, ncol=TP )
        res1 <- sirt_rcpp_gom_em_calcpost( DAT2=dat2, DAT2RESP=dat2.resp, PROBS=probsM,
                                DAT2IND=dat2.ind, PIK=pi.k, KK1=1, weights=weights )
        f.yi.qk <- res1$fyiqk
        f.qk.yi <- res1$f.qk.yi
        pi.k <- res1$pi.k
        n.ik <- array( res1$n.ik, dim=dim(n.ik) )
        N.ik <- res1$N.ik

        ## distribution smoothing for GOMnormal
        if (model=="GOMnormal"){
            res <- stats::cov.wt(x=theta_grid, wt=pi.k, method="ML")
            mu <- res$center
            Sigma <- res$cov
            #- readjust theta_grid
            if (theta_adjust){
                fac1_old <- fac1
                mth <- max(theta0)
                sig_max <- max( sqrt(diag(Sigma)))
                fac1 <- sig_max * ( 6 / mth    )
                theta_grid <- fac1*theta_grid0
                theta.k <- gom_em_normal_to_membership_scores(theta_grid=theta_grid, K=K, TP=TP)
            }
            pi.k <- sirt_dmvnorm_discrete(x=theta_grid, mean=mu, sigma=Sigma)
        }

        # maximize lambda
        if (model %in% c("GOM","GOMnormal") ){
            res <- gom_em_est_lambda( lambda=lambda, I=I, K=K, n.ik=n.ik,
                        numdiff.parm=numdiff.parm, max.increment=max.increment,
                        theta.k=theta.k, msteps=msteps, mstepconv=mstepconv,
                        lambda_partable=lambda_partable, eps=.001, progress=progress )
            lambda <- res$lambda
            se.lambda <- res$se.lambda
            max.increment <- max( abs(lambda-lambda0))/1.2
        }
        if (model=="GOMRasch"){
            res <- gom_em_est_b( lambda=lambda, I=I, K=K, n.ik=n.ik, b=b, theta0.k=theta0.k,
                        numdiff.parm=numdiff.parm, max.increment=max.increment, theta.k=theta.k,
                        msteps=msteps, mstepconv=mstepconv, eps=.001, progress=progress )
            b <- res$b
            se.b <- res$se.b
            lambda <- t( stats::plogis( outer( theta0.k, b, "-" ) ) )
            max.increment <- max( abs(b-b0))/1.2
        }
        utils::flush.console()

        #-- calculate deviance
        ll <- gom_em_compute_total_loglikelihood(f.yi.qk=f.yi.qk, pi.k=pi.k,
                            weights=weights)
        dev <- -2*ll
        # convergence criteria
        conv <- max( abs(lambda-lambda0))
        iter <- iter+1
        devchange <- abs( ( dev - dev0 )/dev0 )

        #**** print progress
        if (progress){
            cat( paste( "   Deviance", "=", round( dev, 4 ),
                if (iter > 1 ){ " | Deviance change=" } else {""},
                if( iter>1){round( - dev + dev0, 6 )} else { ""},"\n",sep="") )
            cat( paste( "    Maximum lambda parameter change=",
                    paste( round(max(abs(lambda-lambda0)),6), collapse=" " ), "\n", sep=""))
            cat( paste( "    Maximum distribution parameter change=",
                    paste( round(max(abs(pi.k-pi.k0)),6), collapse=" " ), "\n", sep=""))
            if (model=="GOMRasch"){
                cat( paste( "    Maximum b parameter change=",
                    paste( round(max(abs(b-b0)),6), collapse=" " ), "\n", sep=""))
            }
        }
    } #--- end EM algorithm

    #--- Newton-Raphson steps
    if (newton_raphson){
        lambda_logit <- as.vector( stats::qlogis(lambda) )
        lambda_logit <- stats::aggregate(lambda_logit, list(lambda_partable$par_index), mean)[,2]
        ind_lambda <- seq(1, max(lambda_partable$par_index))
        names(lambda_logit) <- lambda_partable$par_name[ lambda_partable$free ]
        pi_k_logit <- sirt_probs_to_logit(y=pi.k)
        names(pi_k_logit) <- paste0("pi_Cl", 1:(K-1) )
        ind_pi <- NULL
        ind_mu <- NULL
        ind_sigma <- NULL
        if (model=="GOM"){
            ind_pi <- max(ind_lambda) + seq(1,TP-1)
            x0 <- c(lambda_logit, pi_k_logit)
        }
        if (model=="GOMnormal"){
            K1 <- K-1
            ind_mu <- max(ind_lambda) + seq(1,K1)
            names(mu) <- paste0("mu_",1:K1)
            NC <- K1 * (K1+1) / 2
            ind_sigma <- max(ind_mu) + seq(1,NC)
            Sigma_chol <- t(chol(Sigma))
            x1 <- Sigma_chol[!upper.tri(Sigma_chol)]
            names(x1) <- paste0("sigmachol_",1:NC)
            Sigma1 <- Sigma_chol %*% t(Sigma_chol)
            x0 <- c(lambda_logit, mu, x1)
            names(x0) <- NULL
        }

        #-- optimization
        NP <- length(x0)
        bound <- 8
        lower <- rep(-bound,NP)
        upper <- rep(bound, NP)
        if (model=="GOMnormal"){
            bound2 <- 10
            lower[ c(ind_mu, ind_sigma) ] <- - bound2
            upper[ c(ind_mu, ind_sigma) ] <- bound2
        }

        if (lambda_partable$item_con[1]==0){
            grad <- gom_em_loglike_grad
        } else {
            grad <- NULL
        }
        x0 <- sirt_squeeze(x=x0, lower=lower, upper=upper)

        if (progress){
            cat("\n\n*** Switching to Newton-Raphson Optimization ***\n\n")
        }

        #- optimization
        res <- sirt_optimizer(optimizer=optimizer, par=x0, fn=gom_em_loglike_opt_fun,
                    grad=grad, control=list(trace=2), lower=lower, upper=upper,
                    ind_lambda=ind_lambda, ind_pi=ind_pi, I=I, K=K, theta.k=theta.k,
                    theta0.k=theta0.k, dat2=dat2, dat2.resp=dat2.resp, TP=TP,
                    ind_mu=ind_mu, ind_sigma=ind_sigma, theta_grid=theta_grid,
                    model=model, weights=weights, lambda_partable=lambda_partable,
                    lambda_index_blocks=lambda_index_blocks)
        optimization <- res
        par <- res$par

        #- convert parameters
        res <- gom_em_loglike_parameter_conversion(x=par, ind_lambda=ind_lambda,
                        ind_pi=ind_pi, I=I, K=K, ind_mu=ind_mu, ind_sigma=ind_sigma,
                        theta_grid=theta_grid, model=model, lambda_partable=lambda_partable)
        lambda <- res$lambda
        pi.k <- res$pi.k
        mu <- res$mu
        Sigma <- res$Sigma

        #- recalculate posterior
        res <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=b,
                        theta0.k=theta0.k )
        probs <- res$probs
        probs <- problong2probarray( probres=probs, I=I, TP=TP )
        probsM <- matrix( aperm( probs, c(2,1,3) ), nrow=I*2, ncol=TP )
        res1 <- sirt_rcpp_gom_em_calcpost( DAT2=dat2, DAT2RESP=dat2.resp, PROBS=probsM,
                            DAT2IND=dat2.ind, PIK=pi.k, KK1=1, weights=weights)
        f.yi.qk <- res1$fyiqk
        f.qk.yi <- res1$f.qk.yi
        n.ik <- array( res1$n.ik, dim=dim(n.ik) )
        N.ik <- res1$N.ik
        ll <- gom_em_compute_total_loglikelihood(f.yi.qk=f.yi.qk, pi.k=pi.k,
                            weights=weights)
        dev <- -2*ll
    }

    #--------------- arrange OUTPUT

    #--- Information criteria
    ic <- gom_em_ic( dev=dev, dat2=dat2, I=I, K=K, TP=TP, model=model, weights=weights )

    #--- item parameters
    item <- gom_em_item_parameters( dat2=dat2, dat2.resp=dat2.resp, model=model,
                    b=b, lambda=lambda, K=K, weights=weights, progress=progress )

    EAP.rel <- NULL
    person <- NULL
    if ( model=="GOMRasch"){
        res <- gom_em_est_covariance( f.qk.yi=f.qk.yi, Sigma=Sigma,
                    theta.kM=theta.kM, N=N )
        mu <- res$mu
        Sigma <- res$Sigma
        pi.k <- res$pi.k

        #--- distribution parameters
        c1 <- stats::cov2cor(Sigma)
        if (progress){
            cat("*********************************\n")
            cat("Trait Distribution (Location, Variability)\n")
            cat( " Means: ", round( mu, 3 ), "\n")
            cat( " Standard deviations: ", round( sqrt(diag(Sigma)), 3 ), "\n")
            cat( " Correlation ", round( c1[lower.tri(c1)], 3 ), "\n")
            flush.console()
        }
        # person parameters
        pers <- .smirt.person.parameters( data=dat2, D=2, theta.k=theta.kM,
            p.xi.aj=f.yi.qk, p.aj.xi=f.qk.yi, weights=rep(1,N) )
        person <- pers$person
        EAP.rel <- pers$EAP.rel
        if (progress){
            cat("*********************************\n")
            cat("EAP Reliability=", round(EAP.rel,3), "\n")
        }
    }

    #--- MAP
    MAP <- sirt_MAP(post=f.qk.yi, theta=theta.k)

    #--- EAP
    EAP <- sirt_EAP(post=f.qk.yi, theta=theta.k)
    colnames(EAP) <- paste0("Class",1:K)

    #--- descriptives of classes
    score <- rowSums( weights*dat2 *dat2.resp ) / rowSums( weights*dat2.resp )
    plmat <- FALSE
    if ( is.vector(problevels) ){
        PL <- length(problevels)
    }
    if ( is.matrix(problevels) ){
        PL <- nrow(problevels)
        plmat <- TRUE
    }
    theta.kk0 <- theta.k
    if ( model=="GOMRasch"){
        PL <- 5
        problevels <- seq(0,1,len=PL)
        problevels2 <- c( problevels - diff(problevels)[1]/2, 1.2 )
        for (kk in 1:K){
            for (ll in 1:PL){
                ind.kk <- which( ( theta.k[, kk ] > problevels2[ll] ) &
                    ( theta.k[, kk ] <=problevels2[ll+1] ) )
                theta.kk0[ ind.kk, kk ] <- problevels[ll]
            }
        }
        pi.k <- pi.k[,1]
    }  # end GOMRasch

    classdesc <- NULL
    if ( model=="GOMnormal"){
        plmat <- TRUE
    }
    if ( ! plmat ){
        classdesc <- data.frame( matrix( 0, 2*PL+1, K ) )
        classdesc[1,] <- colSums( theta.kk0 * pi.k )
        for (kk in 1:K){
            for (ll in 1:PL){
                ll1 <- problevels[ll]
                classdesc[ ll+1, kk ] <- sum( pi.k * ( theta.kk0[,kk]==ll1 ) )
            }
        }
        rownames(classdesc)[1] <- "p.Class"
        colnames(classdesc) <- paste0("Class", 1:K )
        rownames(classdesc)[2:(PL+1)] <- paste0( "p.problevel", round( problevels, 3 ), ".class" )
        rownames(classdesc)[PL+2:(PL+1)] <- paste0( "M.problevel", round( problevels, 3 ), ".class" )
        for (kk in 1:K){
            ll.kk <- problevels
            for (ll in 1:PL){
                ind.ll <- which( theta.kk0[,kk]==problevels[ll] )
                rll.wt <- rowSums( f.qk.yi[, ind.ll, drop=FALSE] )
                ll.kk[ll] <- weighted.mean( score, rll.wt )
            }
            classdesc[ PL + 2:(PL+1), kk ]  <- ll.kk
        }
    }

    #--- output
    s2 <- Sys.time()
    res <- list(deviance=dev, ic=ic, item=item, person=person, EAP.rel=EAP.rel,
        MAP=MAP, EAP=EAP, classdesc=classdesc, lambda=lambda, se.lambda=se.lambda,
        mu=mu, Sigma=Sigma, b=b, se.b=se.b, f.yi.qk=f.yi.qk, f.qk.yi=f.qk.yi,
        probs=probs, n.ik=n.ik, iter=iter, dat=dat0, dat2=dat2, dat2.resp=dat2.resp,
        I=I, K=K, TP=TP, G=1, theta.k=theta.k, pi.k=pi.k, problevels=problevels,
        model=model, plmat=plmat, mu=mu, Sigma=Sigma, weights=weights,
        lambda.index=lambda.index, lambda_partable=lambda_partable,
        optimization=optimization, newton_raphson=newton_raphson,
        s1=s1, s2=s2, time_diff=s2-s1, CALL=CALL)
    class(res) <- "gom"
    return(res)
}


# a1 <- Sys.time(); cat("probs \n"); print(a1-a0); a0 <- a1
