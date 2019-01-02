## File Name: noharm.sirt.R
## File Version: 0.863


########################################
# NOHARM implementation in R
########################################

#---------------
# NOHARM
# Normal-Ogive Harmonic Analysis Robust Method
#---------------
# The other function with linear and loglinear models of parameters
# is denoted as 'noharm.loglm.sirt'.
#---------------

noharm.sirt <- function(dat, weights=NULL, Fval=NULL, Fpatt=NULL,
    Pval=NULL, Ppatt=NULL, Psival=NULL, Psipatt=NULL, dimensions=NULL,
    lower=rep(0,ncol(dat)), upper=rep(1,ncol(dat)), wgtm=NULL,
    modesttype=1, pos.loading=FALSE, pos.variance=FALSE,
    pos.residcorr=FALSE, maxiter=1000, conv=1e-6,  increment.factor=1.01,
    reliability=TRUE )
{
    s1 <- Sys.time()
    CALL <- match.call()
    #*** data preprocessing
    res <- noharm_sirt_preproc( dat=dat, weights=weights, Fpatt=Fpatt, Fval=Fval,
                Ppatt=Ppatt, Pval=Pval, Psipatt=Psipatt, Psival=Psival, wgtm=wgtm,
                dimensions=dimensions )
    ss <- res$ss
    pm0 <- res$pm0
    sumwgtm <- res$sumwgtm
    model.type <- res$model.type
    modtype <- res$modtype
    N <- res$N
    I <- res$I
    dat0 <- res$dat0
    Fpatt <- res$Fpatt
    Fval <- res$Fval
    Ppatt <- res$Ppatt
    Pval <- res$Pval
    Psipatt <- res$Psipatt
    Psival <- res$Psival
    D <- res$D
    pm <- res$pm
    dat <- res$dat
    dat.resp <- res$dat.resp
    wgtm <- res$wgtm
    F_dimnames <- res$F_dimnames
    items <- res$items
    wgtm.default <- res$wgtm.default

    #    modesttype <- 1        # NOHARM estimation
    #    modesttype <- 2        # estimation using tetrachoric correlations
    if (modesttype==2){
        pm2 <- tetrachoric2( dat, method="Bo")
        pm <- pm2$rho
        betaj <- pm2$tau
        lower <- 0*lower
        upper <- 1+lower
    }
    #----
    # compute betaj
    if (modesttype==1){
        betaj <- - stats::qnorm( ( diag(pm) - lower ) / ( upper - lower ) )
    }
    # compute bjk coefficients
    # include lower and upper asymptotes here
    b0 <- lower + (upper-lower) * stats::pnorm(-betaj)
    b1 <- (upper-lower) * stats::dnorm(betaj)
    b2 <- betaj * b1 / sqrt(2)
    b3 <- ( betaj^2 -1 ) * b1 / sqrt(6)

    # create fixed cofficients
    b0.jk <- noharm_sirt_outer_coefs(x=b0)
    b1.jk <- noharm_sirt_outer_coefs(x=b1)
    b2.jk <- noharm_sirt_outer_coefs(x=b2)
    b3.jk <- noharm_sirt_outer_coefs(x=b3)

    parchange <- 1
    changeF <- changeP <- changePsi <- 0
    iter <- 0
    maxincrement <- .2
    estPsi <-  sum( Psipatt > 0 ) > 0
    estF <- sum( Fpatt > 0 ) > 0
    estP <- sum( Ppatt > 0 ) > 0
    estpars <- list(estF=estF, estP=estP, estPsi=estPsi)
    eps <- 2*conv

    #- list of arguments in optimization routine
    args <- list( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I, D=D,
                        b0jk=b0.jk, b1jk=b1.jk, b2jk=b2.jk, b3jk=b3.jk, wgtm=wgtm, pm=pm, Psival=Psival,
                        Psipatt=Psipatt, maxincrement=maxincrement, modtype=modesttype )

    #**** begin algorithm
    while( ( iter < maxiter ) & ( parchange > conv ) ){
        maxincrement <- maxincrement  / increment.factor
        args$maxincrement <- maxincrement

        #---- update F
        if (estF){
            res <- do.call(what=sirt_rcpp_noharm_estF, args=args)
            changeF <- res$change
            Fval <- res$Fval_
            # F_der <- res$F_der
            if ( pos.loading ){
                Fval[ Fval < 0 ] <- eps
            }
            args$Fval <- Fval
        }
        #---- update P
        if (estP){
            res <- do.call(what=sirt_rcpp_noharm_estP, args=args)
            changeP <- res$change
            Pval <- res$Pval_
            if ( pos.variance ){
                diag(Pval)[ diag(Pval) < 0 ] <- eps
            }
            diag(Pval)[ is.na(diag(Pval)) ] <- eps
            args$Pval <- Pval
        }
        #---- update Psi
        if (estPsi){
            res <- do.call(what=sirt_rcpp_noharm_estPsi, args=args)
            changePsi <- res$change
            Psival <- res$Psival_
            if ( pos.residcorr ){
                Psival[ Psival < 0 ] <- eps
            }
            args$Psival <- Psival
        }
        parchange <- max( c(changeP, changeF, changePsi) )
        iter <- iter + 1
    }
    #**** end algorithm

    #--- include colnames in estimated matrices
    colnames(Fval) <- F_dimnames
    rownames(Fval) <- items
    rownames(Pval) <- colnames(Pval) <- F_dimnames
    rownames(Psival) <- colnames(Psival) <- items

    #* compute final constants
    res <- noharm_sirt_compute_final_constants( Fval=Fval, Pval=Pval, betaj=betaj,
                modesttype=modesttype )
    f0 <- res$f0
    uqn <- res$uqn
    loadingsF <- res$loadingsF

    #* compute residuals
    res <- noharm_sirt_est_residuals( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I,
                D=D, b0.jk=b0.jk, b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk, wgtm=wgtm, pm=pm,
                Psival=Psival, Psipatt=Psipatt, modesttype=modesttype )
    residuals <- res$residuals
    rmsr <- res$rmsr
    tanaka <- res$tanaka

    #**** arrange output list
    res <- list( tanaka=tanaka, rmsr=rmsr, N.itempair=ss, pm=pm0, wgtm=wgtm, sumwgtm=sumwgtm,
                lower=lower, upper=upper, residuals=residuals, final.constants=f0,
                factor.cor=Pval, thresholds=betaj, uniquenesses=uqn, loadings=loadingsF,
                loadings.theta=Fval, residcorr=Psival, model.type=model.type, modtype=modtype,
                Nobs=N, Nitems=I, Fpatt=Fpatt, Ppatt=Ppatt, Psipatt=Psipatt,
                dat=dat0, systime=Sys.time(), dimensions=D, display.fit=5, iter=iter )

    #*** number of estimated parameters
    Nestpars <- noharm_sirt_number_estimated_parameters( I=I, Fpatt=Fpatt,
                    Ppatt=Ppatt, Psipatt=Psipatt )
    res$Nestpars <- Nestpars

    #*** compute chi square statistics
    res <- noharm_sirt_compute_chi_square_statistics( res=res,
                residuals=residuals, pm=pm, N=N, I=I, sumwgtm=sumwgtm, modesttype=modesttype,
                Nestpars=Nestpars )

    #*** compute Green-Yang reliability
    v1 <- 1 * ( ( wgtm - diag(wgtm) ) > 0 )
    # ensure positive definiteness
    res$omega.rel <- NA
    L0 <- sqrt( diag(Pval) )
    N1 <- length(L0)
    L1inv <- L1 <- matrix( 0, nrow=N1, ncol=N1 )
    rownames(L1) <- rownames(L1inv) <- colnames(L1) <- colnames(L1inv) <- rownames(Pval)
    diag(L1) <- L0
    diag(L1inv) <- 1 / L0
    Fval2 <- Fval %*% L1
    Pval2 <- L1inv %*% Pval %*% L1inv
    dj <- sqrt( diag( Fval2 %*% Pval2 %*% t(Fval2) ) )
    ej <- sqrt( 1+dj^2 )
    Fval2 <- Fval2 / ej
    standardized.solution <- list( Fval=Fval2, Pval=Pval2 )
    if (reliability){
        res$omega.rel <- reliability.nonlinearSEM(facloadings=Fval2,
            thresh=res$thresholds, resid.cov=res$residcorr, cor.factors=Pval2 )$omega.rel
    }
    if ( sum(v1) + I < I^2 ){
        res$omega.rel <- NA
    }

    #*** rotated solution for EFA
    if (model.type=="EFA"){
        res <- noharm_sirt_efa_rotated_solution( res=res, items=items,
                    F_dimnames=F_dimnames )
    }
    if (modesttype==2){
        res$tetracor <- pm
    }

    #***----- more output
    res$estpars <- estpars
    res$modesttype <- modesttype
    res$guesses <- res$lower
    res$wgtm.default <- wgtm.default
    res$standardized.solution <- standardized.solution
    s2 <- Sys.time()
    res$s1 <- s1
    res$s2 <- s2
    res$CALL <- CALL
    class(res) <- 'noharm.sirt'
    return(res)
}
