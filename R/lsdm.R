## File Name: lsdm.R
## File Version: 1.387


# LSDM - Least Squares Distance Method
# LSDM -- Least Squares Distance Method of Cognitive Validation
# Reference: Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.
lsdm <- function( data, Qmatrix, theta=seq(-3,3,by=.5), wgt_theta=rep(1, length(theta)),
            distance="L2", quant.list=c(.5,.65,.8), b=NULL, a=rep(1,nrow(Qmatrix)),
            c=rep(0,nrow(Qmatrix)) )
{
    CALL <- match.call()
    s1 <- Sys.time()
    # generate sequence for display
    display.separate <- paste( rep(".", each=65 ), collapse="" )
    # display progress
    cat(display.separate, "\n" )
    cat( "LSDM -- Least Squares Distance Method of Cognitive Validation \n")
    cat("Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.\n")
    cat(display.separate, "\n")
    TP <- length(theta)
    I <- nrow(Qmatrix)
    thetaM <- sirt_matrix2(theta, nrow=I)
    wgt_theta <- sirt_sum_norm(x=wgt_theta)
    if (! is.null(b) ){
        eins <- rep(1, TP)
        cM <- matrix(c, nrow=I, ncol=TP)
        bM <- matrix(b, nrow=I, ncol=TP)
        aM <- matrix(a, nrow=I, ncol=TP)
        data <- cM + ( 1 - cM )* stats::plogis( aM * ( thetaM - bM ) )
    }
    Qmatrix <- as.matrix(Qmatrix)

    # print Q matrix
    cat("\nQmatrix\n\n")
    cmax <- apply( Qmatrix, 2, max )
    Qmatrix <- Qmatrix / sirt_matrix2(cmax, nrow=I)
    print(Qmatrix)
    cat("\n")
    d1 <- det( crossprod(Qmatrix) )
    # warning for singular Q matrices
    if (abs(d1) < 1E-8){
        stop("You used a singular Q matrix as an input. LSDM cannot be computed.\n")
    }
    est.icc <- TRUE
    I <- nrow(data)
    L <- ncol(data)
    K <- ncol(Qmatrix)
    if ( is.null(rownames(Qmatrix)) ){
        rownames(Qmatrix) <- rownames(data)
    }
    # log probability functions
    data1 <- data
    eps <- 1e-3
    logdata <- log( data1 + eps)
    # estimate item parameter and item quantiles
    cat("Estimation of Item Parameters \n")
    utils::flush.console()
    icc.pars <- lsdm_est_logist_quant( probcurves=data, theta=theta,
                    quantiles=quant.list, wgt_theta=wgt_theta, est.icc=est.icc,
                    b=b, a=a)
    cat(display.separate, "\n" )

    #***** Estimate attribute response curves
    cat("Estimation of Attribute Parameters \n")
    utils::flush.console()

    Q <- Qmatrix
    eps_obj <- 1e-4
    if (distance=="L2"){ pow_obj <- 2 }
    if (distance=="L1"){ pow_obj <- 1 }

    #** define optimization function for LSDM method
        # P=A1 * A2
        # log(P)=log(A1) + log(A2) for every tt
        #- apply restricted linear model
    lsdm_obj_fn <- function(x, y_tt){
        dist_log <- y_tt - Q %*% x
        if (pow_obj==2){
            val <- sum(dist_log^2)
        }
        if (pow_obj==1){
            val <- sum(sirt_abs_smooth(x=dist_log, eps=eps_obj))
        }
        return(val)
    }
    lsdm_obj_gr <- function(x, y_tt){
        K <- length(x)
        I <- length(y_tt)
        dist_log <- y_tt - Q %*% x
        dist_log <- matrix(dist_log, nrow=I, ncol=K)
        if (pow_obj==2){
            grad <- -2*colSums(dist_log*Q)
        }
        if (pow_obj==1){
            der1 <- sirt_abs_smooth(x=dist_log, eps=eps_obj, deriv=1 )
            grad <- -colSums( der1*Q )
        }
        return(grad)
    }

    log.arc0 <- matrix(NA, nrow=K, ncol=L)
    par0 <- rep(-2,K)
    upper <- rep(0,K)
    for (tt in 1:L){
        y_tt <- logdata[,tt]
        res <- stats::optim(par=par0, fn=lsdm_obj_fn, gr=lsdm_obj_gr, method="L-BFGS-B",
                        upper=upper, y_tt=y_tt)
        par0 <- res$par
        log.arc0[,tt] <- par0
    }

    # exponentiate attribute response curve
    arc0 <- exp(log.arc0)
    rownames(arc0) <- colnames(Qmatrix)

    #--- estimate "ordinary" LLTM
    lltm.res1 <- stats::lm( as.numeric(icc.pars$b.1PL) ~ 0 + as.matrix(Qmatrix) )
    slltm.res1 <- summary(lltm.res1)
    cat(display.separate, "\n")

    # calculate Rasch data predicted by LLTM
    b_lltm <- lltm.res1$fitted
    data.lltm <- stats::plogis(thetaM - b_lltm)
    # estimate attribute parameter and attribute quantiles
    arc0.pars <- lsdm_est_logist_quant( probcurves=arc0, theta=theta,
                quantiles=quant.list, est.icc=est.icc, wgt_theta=wgt_theta )
    arc0.pars$eta.LLTM <- coef(lltm.res1)
    arc0.pars$se.LLTM <- slltm.res1[[4]][,2]
    arc0.pars$pval.LLTM <- slltm.res1[[4]][,4]
    W <- matrix( NA, nrow=I, ncol=K )
    for (ii in 1:I){
        index.ii <- which( Qmatrix[ii,] > 0 )
        L.ii <- length(index.ii)
        x.ii <- t( log.arc0[ index.ii, ] )
        x.ii <- matrix( x.ii, ncol=length(index.ii) )
        y.ii <- as.numeric(logdata[ii,])
        mod1.ii <- stats::lm( y.ii ~ 0 + x.ii )
        W[ii,index.ii] <- coef( mod1.ii )
    }

    #---- evaluate goodness of fit
    data0.fitted <- exp( Qmatrix %*% log.arc0 )
    rownames(data.lltm) <- rownames(data0.fitted) <- rownames(data)
    # MAD for original model (Dimitrov)
    mad0 <- lsdm_irf_distance_mad(data=data, data_fitted=data0.fitted,
                            wgt_theta=wgt_theta, use_abs=TRUE)
    md0 <- lsdm_irf_distance_mad(data=data, data_fitted=data0.fitted,
                            wgt_theta=wgt_theta, use_abs=FALSE)
    mm0 <- mean(mad0)

    mad.lltm <- lsdm_irf_distance_mad(data=data, data_fitted=data.lltm,
                            wgt_theta=wgt_theta, use_abs=TRUE)
    md.lltm <- lsdm_irf_distance_mad(data=data, data_fitted=data.lltm,
                            wgt_theta=wgt_theta, use_abs=FALSE)
    mm.lltm <- mean(mad.lltm)

    #* model fit LSDM
    cat(paste( "Model Fit LSDM   -  Mean MAD:",
                    formatC( round( mm0, 3 ),digits=3, width=6),
                    "    Median MAD:", formatC( round( median(mad0), 3 ),
                                digits=3, width=6), "\n") )
    cat(paste( "Model Fit LLTM   -  Mean MAD:",
                    formatC( round( mm.lltm, 3 ),digits=3, width=6),
                    "    Median MAD:", formatC( round( median(mad.lltm), 3 ),
                                    digits=3, width=6),
                    "   R^2=", format( round( slltm.res1$r.squared, 3 ),digits=3), "\n") )
    item <- data.frame( N.skills=rowSums(Qmatrix), mad.lsdm=mad0,
                        mad.lltm=mad.lltm, md.lsdm=md0, md.lltm=md.lltm, icc.pars )
    colnames(W) <- colnames(Qmatrix)
    rownames(W) <- rownames(Qmatrix)

    #--- output
    s2 <- Sys.time()
    time <- list(s1=s1, s2=s2, time_diff=s2-s1)
    res <- list( mean.mad.lsdm0=mm0,  mean.mad.lltm=mm.lltm, attr.curves=arc0,
                attr.pars=arc0.pars, data.fitted=data0.fitted, theta=theta,
                item=item, data=data, Qmatrix=Qmatrix, lltm=lltm.res1, W=W,
                distance=distance, CALL=CALL, time=time )
    class(res) <- "lsdm"
    return(res)
}

