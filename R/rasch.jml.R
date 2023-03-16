## File Name: rasch.jml.R
## File Version: 3.315


rasch.jml <- function( dat, method="MLE", b.init=NULL, constraints=NULL, weights=NULL,
                    center="persons", glob.conv=10^(-6), conv1=0.00001, conv2=0.001, progress=TRUE,
                    bsteps=4, thetasteps=2, wle.adj=0, jmliter=100, prox=TRUE, proxiter=30, proxconv=0.01,
                    dp=NULL, theta.init=NULL, calc.fit=TRUE, prior_sd=NULL)
{
    CALL <- match.call()
    s1 <- Sys.time()
    # warning messages
    if ( max( dat, na.rm=TRUE )  > 1){
        stop( "You should only have dichotomous data")
    }
    # zz0 <- Sys.time()

    # display
    if ( progress ){
        cat("-----------------------------------------------------------------\n")
        cat("Joint Maximum Likelihood Estimation \n")
        cat("Rasch Model \n")
        cat("-----------------------------------------------------------------\n")
    }

    # prior for ability
    if (! is.null(prior_sd) ){
        is_prior <- TRUE
    } else {
        is_prior <- FALSE
    }

    #--- centering of parameters
    centeritems <- FALSE
    centerpersons <- FALSE
    if (center=="persons"){
        centerpersons <- TRUE
    }
    if (center=="items"){
        centeritems <- TRUE
    }
    if (! is.null(constraints) ){
        prox <- FALSE
        centerpersons <- FALSE
        center <- "none"
    }

    # data preparations
    use_weights <- FALSE
    if (!is.null(weights)){
        use_weights <- TRUE
        prox <- FALSE
    }
    if ( is.null(dp) ){
        dp <- data.prep( dat, weights=weights, use.freqpatt=FALSE,
                standardize_weights=FALSE)
    }
    dat1 <- dp$dat1
    dat2 <- dp$dat2
    dat2.resp <- dp$dat2.resp
    freq.patt <- dp$freq.patt
    I <- dp$I
    N <- n <- dp$n

    # exclude extreme item response pattern (if method==MLE)
    if ( ( method=="MLE" ) & ( ! is_prior ) ){
        ind <- which( ! dat1$mean %in% c(0,1) )
        dat1 <- dat1[ ind,]
        dat2 <- dat2[ ind,]
        dat2.resp <- dat2.resp[ ind,]
        N <- nrow(dat1)
    }

    #**** exist some missings?
    some.missings <- 1 * ( mean( rowMeans( dat2.resp) ) < 1 )
    # compute missing response pattern
    dat1$mp <- paste("MP", paste(rep(1,I), collapse=""), sep="")
    dat1$mp.index <- 1
    if ( some.missings ){
        dat1$mp.index <- resp.pattern2( dat2.resp )$mp.index
    }
    dat1$sumscore <- rowSums( dat2 * dat2.resp )
    ti <- paste( dat1$mp.index, "-", dat1$sumscore, sep="")
    dat1$theta.index <- match( ti, unique( ti ))
    dat1$caseid <- seq( 1, nrow(dat1) )

    theta.pattern <- as.data.frame( stats::aggregate( dat1$caseid, list(dat1$theta.index), min ))
    colnames(theta.pattern) <- c("theta.index", "caseid")

    theta.pattern[,"Freq"] <- rowsum( dat1$Freq, dat1$theta.index )[,1]
    TPsapp <- 1:nrow(theta.pattern)
    freq.thetapattern <- rowsum( dat1$Freq, dat1$theta.index )[,1]
    if ( some.missings ){
        freq.dat.resp.thetapattern <- rowsum( dat1$Freq  * dat2.resp, dat1$theta.index )
    } else {
        freq.dat.resp.thetapattern <- matrix( freq.thetapattern, nrow=length(freq.thetapattern), ncol=I )
    }

    #** compute sufficient statistic
    suffB <- colSums(  - dat1[,2] * dat2.resp * dat2 )
    d2d2r <- dat2 * dat2.resp
    d2d2m <- ( 1 - dat2 ) * dat2.resp

    # starting values item parameter
    if ( is.null(b.init) ){
        b <- - stats::qlogis( colSums( dat2 * dat2.resp * dat1[,2] ) / colSums( dat1[,2] * dat2.resp ) )
    } else {
        b <- b.init
    }
    if (!is.null(constraints)){
        b[ constraints[,1] ] <- constraints[,2]
    }

    # starting values person paramters
    if ( is.null(theta.init)){
        theta <- stats::qlogis( dat1$mean * ( 1  - 1 / I ) + 1 / ( 2*I ) )
    } else {
        theta <- theta.init
        theta[ theta==Inf ] <- 20
        theta[ theta==-Inf ] <- -20
    }

    #--- PROX estimation
    if (prox){
        if (method=="WLE" ){
            ind <- which( ! dat1$mean %in% c(0,1) )
        } else {
            ind <- 1:( nrow(dat1 ))
        }
        prox.res <- rasch.prox( dat=dat2[ind,], dat.resp=dat2.resp[ind,], freq=dat1[ind,2],
                        conv=proxconv, progress=progress, maxiter=proxiter )
        b <- prox.res$b
        theta[ind] <- prox.res$theta
    }
    # centering
    theta <- rasch_jml_centerpersons(theta=theta, dat1=dat1, centerpersons=centerpersons)
    b <- rasch_jml_centeritems(b=b, centeritems=centeritems)

    # initial settings
    dev0 <- 1
    dev.change <- 1
    iter <- 0
    b0 <- 0
    bconv <- 10^10
    disp <- "...............................\n"

    #------------------------------------------------------
    # JML Iteration Algorithm
    while ( ( ( dev.change > glob.conv ) & ( iter < jmliter ) ) | ( bconv > conv1 ) ){
        b0 <- b
        # update item difficulties
        if ( progress ){
            cat(disp)
            cat("JML Iteration", iter +1, "\n")
            cat("  Item parameters |")
        }
        b <- rasch_jml_update_b( b0, theta=theta[theta.pattern$caseid],
                        freq.thetapattern=freq.thetapattern,
                        freq.dat.resp.thetapattern=freq.dat.resp.thetapattern,
                        constraints=constraints, conv=conv2,
                        suffB=suffB, progress=progress, bsteps=bsteps)
        b <- rasch_jml_centeritems(b=b, centeritems=centeritems)
        # ability estimate (method==MLE or method==WLE)
        ind <- theta.pattern$caseid
        if (method=="WLE" ){
            m1 <- wle.rasch( dat=dat2[ind,], b=b, theta=theta[ind],
                        dat.resp=dat2.resp[ind, ], conv=conv2, progress=progress,
                        maxit=thetasteps, wle.adj=wle.adj )
        }
        if (method=="MLE" ){
            m1 <- mle.rasch( dat=dat2[ind,], b=b, theta=theta[ind], dat.resp=dat2.resp[ind,],
                            conv=conv2,    progress=progress, prior_sd=prior_sd )
        }
        theta <- m1$theta
        theta <- theta[ dat1$theta.index ]
        # centering
        theta <- rasch_jml_centerpersons(theta=theta, dat1=dat1, centerpersons=centerpersons)
        # calculate Log-Likelihood and Deviance
        p.ia <- (m1$p.ia)[ dat1$theta.index, ]
        dev1 <- -2*sum( rowSums( ( d2d2r*log(p.ia) + d2d2m*log( 1 - p.ia ))) * dat1[,2] )
        # relative deviance change
        dev1a <- - (dev1 - dev0)
        dev.change <- abs( dev1 - dev0 ) / dev0 ;
        dev0 <- dev1
        # iteration index
        iter <- iter + 1
        # display convergence
        bconv <- max( abs(b - b0 ) )
        if ( progress ){
            cat( paste( "\n  Deviance", sirt_sign_space(), round( dev1, 5 ) ))
            if (iter > 1){
                cat( paste0( " | Deviance change", sirt_sign_space(), round( dev1a, 5 ) ) )
            }
            cat( paste0("\n  Max. parm. change", sirt_sign_space(), round( bconv, 6 ), "\n"  ) )
            utils::flush.console()
        }
    }
    #------------- end algorithm

    if ( progress ){
        if ( iter < jmliter ){
            cat( paste( "Convergence reached in", iter, "JML Iterations \n" ) )
        } else {
            cat( paste( "Analysis terminated at", iter, "JML Iterations \n" ) )
        }
        cat("-----------------------------------------------------------------\n\n")
    }
    dat1 <- dat1[,1:3]
    # arrange ability estimates
    abil <- data.frame( dat1, theta )
    abil <- merge( freq.patt, abil[ c(1,4)], 1, 1, all=T )
    abil <- abil[ order(abil[,3]), c( 3,1,2,4) ]
    colnames(abil) <- c("person", "pattern", "meancorrect", "theta" )
    if (method=="MLE"){
        abil[abil$meancorrect==0, 4 ] <- -Inf
        abil[abil$meancorrect==1, 4 ] <- Inf
    }

    # standard error of person parameter
    abil$se.theta <- rasch.info.mle( dat=dat, theta=abil$theta, b=b * (I-1)/I)

    #- summary ability parameters
    theta_summary <- rasch_jml_person_parameters_summary(x=abil[ is.finite(abil$theta), "theta"])

    # empirical discrimination
    a.i <- rasch_jml_emp_discrim( theta=theta, b=b, dat=dat2, dat.resp=dat2.resp, freq=dat1[,2] )

    # calculate itemfit statistics
    if (calc.fit){
        fit <- rasch.itemfit( theta0=abil$theta, b=b, dat )
    } else {
        fit <- NULL
    }


    # standard error for item parameter
    se.b <- .se.item.rasch( theta=theta, b=b, dat.resp=dat2.resp, freq=dat1[,2],
                    constraints=constraints)

    # Data frame for "item side"
    I <- length(b)
    # item parameter constraints and correction formula
    UJJ <- 0
    # item summary
    dfr <- data.frame( N=colSums(dat2.resp*dat1[,2]),
                p=colMeans( dat, na.rm=TRUE ), itemdiff=b,
                itemdiff.correction=b*(I-1)/I, se=se.b, discr=a.i )
    if ( ! is.null(constraints)    ){
        dfr$itemdiff.correction <- dfr$itemdiff
    }


    # include item fit statistics                              #
    if (calc.fit){
        dfr <- data.frame( dfr, fit )
    }
    rownames(dfr) <- colnames(dat)

    # processed data
    data.proc <- list( "data"=dat2, "data.resp"=dat2.resp, "theta"=theta )
    s2 <- Sys.time()
    time_diff <- s2-s1
    #--- output
    res <- list( item=dfr, person=abil, theta_summary=theta_summary,
                    method=method,dat=dat, deviance=dev1,
                    data.proc=data.proc, dp=dp, constraints=constraints, N=N, I=I,
                    method=method, center=center, iter=iter,
                    CALL=CALL, s1=s1, s2=s2, time_diff=time_diff )
    class(res) <- "rasch.jml"
    return(res)
}

# cat("fit") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1

#---------------------------------------------------------------------------------
# Functions for calculating standard errors in the Rasch Model                   #
.se.item.rasch <- function( theta, b, dat.resp, freq, constraints ){
        # standard error for item parameter
        pij <- .prob.rasch( theta=theta, b=b )
        se.b <- sqrt( 1/colSums( freq * pij * ( 1 - pij ) * dat.resp )) # calculate information function
        if (!is.null(constraints)){ se.b[ constraints[,1]] <- NA }
        se.b
        }
#---------------------------------------------------------------------------------

#-----------------------------------------------------
# Calculate of standard measurement error using
# test information function from the Rasch model
##NS export(rasch.info.mle)
rasch.info.mle <- function( dat, theta, b){
    #......................
    # INPUT:
    # dat   ... original data
    # theta ... ability estimates
    # b     ... item difficulties
    dat.resp <- 1 - is.na(dat)
    p.ia <- stats::plogis( outer( theta, b, "-" ) )
    # calculate information function
    info <- rowSums( dat.resp * p.ia * ( 1- p.ia ) )
    # calculate standard error
    sqrt( 1 / info )
    }
#------------------------------------------------------

# auxiliary functions
#-----------------------------------------------------------------------------------------------------------#
# gradient Rasch function
.b.gradient.rasch <- function( b, theta, freq, dat, dat.resp ){
    colSums(  - freq * dat.resp * ( dat - stats::plogis(  outer( theta, b, "-" ) ) ) )
    }
# Information Matrix item difficulties
.infb.rasch1 <- function( b, theta, freq ){
#    p.ia <- plogis( outer( theta, b, "-" ) )
    p.ia <- stats::plogis( theta, matrix( b, nrow=length(theta), length(b), byrow=T ) )
    colSums( - freq * p.ia * ( 1- p.ia ) )
    }
.infb.rasch <- function( b, theta, freq ){
    p.ia <- stats::plogis( outer( theta, b, "-" ) )
    colSums( - freq * p.ia * ( 1- p.ia ) )
    }

#---- calculate P_i ( theta)
.prob.rasch <- function( theta, b )
{
    bM <- matrix( b, nrow=length(theta), length(b), byrow=TRUE )
    res <- stats::plogis( theta - bM )
    return(res)
}

# calculate P_i(theta) for 3PL model
.prob.3pl <- function( theta, b, a, c){
    l1  <- rep( 1, length(theta) )
    cM <- outer( l1, c )
    aM <- outer( l1, a )
    bM <- outer( l1, b )
    thetaM <- outer( theta, rep(1,length(b) ) )
    cM + (1-cM) * stats::plogis( aM*(thetaM - bM ) )
    }
#-----------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------
# calculate information function in the Rasch model
.info.rasch <- function( theta, b ){
    pij <-  .prob.rasch( theta=theta, b=b )
    sqrt(  1 / rowSums( pij * ( 1 - pij ) ) )
        }
#-------------------------------------------------------------------#
# calculate (mean) reliability for known variance and mean for      #
# latent trait distribution                                         #
##NS # export(mle.reliability.rasch)
mle.reliability.rasch <- function(  b, mean.abil=NULL, var.abil=NULL, npoints=200,
            theta=NULL){
    # INPUT:
    # b ... item difficulties
    # mean.abil ... mean of ability distribution
    # var.abil ... variance of ability distribution
    # npoints ... number of discretization points
    # theta grid ... optional
    #...................................
    q.grid <- seq( 1 /npoints, 1 - 1/npoints, len=npoints )
    if ( is.null(theta)){
            theta <- stats::qnorm( q.grid, mean=mean.abil, sd=sqrt( var.abil ) )
                        }
    info <- .info.rasch( theta=theta, b=b )
    # mean estimate of variance
    mean.var.est <- mean( info^2 )
    # Note: info.grid is the standard error
    list( "theta.grid"=theta, "info.grid"=info, "mean.var.est"=mean.var.est,
            "I"=length(b),     "mle.rel"=var.abil /( var.abil + mean.var.est ) )
    }
#-------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------#
# update item difficulty estimation (Rasch model)
.update.b.rasch.jml1 <- function( b, theta, freq, dat, dat.resp,  constraints=NULL, conv=.0001,
                suffB, thetaindex=NULL, progress=FALSE ){
    #-------------------------------------------------------
    # INPUT:
    # b     ... initial item difficulties
    # theta ... ability estimate
    # freq  ... frequency of item response pattern (theta estimate)
    # dat   ... data frame
    # update item difficulties
    b.change <- 1
    while( max( abs( b.change  ) ) > conv ){
        p.ia <- stats::plogis( theta, matrix( b, nrow=length(theta), length(b), byrow=T ) )
        if ( ! is.null( thetaindex ) ){
            p.ia <- p.ia[ thetaindex, ]
                                    }
        deriv <- colSums( - freq * p.ia * ( 1- p.ia ) )
        diff <- suffB + colSums( freq * dat.resp *  p.ia  )
        b.change <-  diff / deriv
        if (! is.null(constraints)){
                b.change[ constraints[,1] ] <- 0
                    }
        if (progress){    cat("-"); utils::flush.console() }
        b <- b - b.change
        }
    b
    }
#-----------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------#
# update item difficulty estimation (Rasch model)
.update.b.rasch.jml <- function( b, theta, freq, dat, dat.resp,  constraints=NULL, conv=.0001 ){
    #-------------------------------------------------------
    # INPUT:
    # b     ... initial item difficulties
    # theta ... ability estimate
    # freq  ... frequency of item response pattern (theta estimate)
    # dat   ... data frame
    # update item difficulties
    b.change <- 1
    while( max( abs( b.change  ) ) > conv ){
        deriv <- .infb.rasch( b, theta, freq )
        b.change <-  .b.gradient.rasch( b, theta, freq, dat, dat.resp ) / deriv
        if (! is.null(constraints)){  b.change[ constraints[,1] ] <- 0 }
        b <- b - b.change
        }
    b
    }
#-----------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------
# calculate empirical item discrimination
item.discrim <- function( dat, score ){
    # INPUT:
    # dat       ... data frame
    # score     ... score with which items are to be correlated (e.g. sum score, WLE, ...)
    apply( dat, 2, FUN=function(variable){ cor( variable,  score, use="complete.obs") } )
    }
#------------------------------------------------------------------------------------------------------




##########################################################################################
rasch.jml.jackknife1 <- function( jmlobj  ){
    #, jackunits=NULL
    jackunits <- NULL
    mod <- jmlobj
    # define jackknife units
    cat("Joint Maximum Likelihood Estimation \nJackknife Estimation\n")
    dat <- jmlobj$dat
    Nobs <- nrow( dat )
    I <- ncol(dat)
    if ( is.null(jackunits) ){ jackunits <- 1:I }
    UJ <- max( jackunits)
    # define jackknife units
    if ( length(jackunits) > 1){
            jackunits <- match( jackunits, unique( jackunits) )
                                }
    # define progress bar
    u.jackunits <- unique( jackunits )
    cat( UJ, "Jackknife Units are used\n" )
    prbar <- paste( "|", paste(rep("-", 20 ), collapse=""), "PROGRESS",
                    paste(rep("-",20 ), collapse=""), "|", sep="")
    cat(prbar,"\n")
    #***************
    # initialize matrices
    cat("|")
    j.itemdiff <- matrix(NA,, nrow=I, ncol=UJ )
    # calculate jackunits indices when progress should be displayed
    prbar.display <- floor( seq( 1, max(u.jackunits), len=48 ) )
    # collect constraints
    constraints <- jmlobj$constraints
    #**************************
    # BEGIN JACKKNIFE
    for ( uu in u.jackunits ){
        # uu <- 10
        dat.uu <- dat[, - uu]
        if ( is.null(constraints) ){
            constraints.uu <- NULL } else {
            constraints.uu <- constraints[ constraints[,1] !=uu, ]
            M1 <- cbind( setdiff( 1:I, uu ), seq(1,ncol(dat.uu) ) )
            constraints.uu <- merge( x=constraints.uu, y=M1, by=1, all.x=TRUE )
            constraints.uu <- constraints.uu[, c(3,2) ]
                        }

        mod.uu <- rasch.jml( dat=dat.uu, b.init=(mod$item$itemdiff)[-uu],
                            prox=FALSE, progress=FALSE,
#                            dp=mod$dp,
                            constraints=constraints.uu,
                            method=mod$method, calc.fit=FALSE )
        j.itemdiff[-uu,uu] <- mod.uu$item$itemdiff
        if ( uu %in% prbar.display ){
                    suu <- sum( prbar.display==uu )
                    cat(paste( rep("-",suu), collapse="")) ; utils::flush.console()
                                }
                    }
    cat("|\n\n")
    # END JACKKNIFE
    #*************************
    dfr <- data.frame( mod$item[,1:2]   )
    # item difficulty
    dfr$b.JML <- mod$item$itemdiff
    dfr$b.JMLcorr <- mod$item$itemdiff.correction
    mean.jack <- rowMeans( j.itemdiff, na.rm=T)
    # calculate correction factor
    UJJ <- UJ
    if ( ! is.null(constraints) ){
            UJJ <- UJJ - nrow(constraints)
                }
    dfr$b.jack <- UJJ * mod$item$itemdiff - (UJJ-1)*mean.jack
    dfr$est.bias <- - dfr$b.jack + dfr$b.JML
    dfr$b.JMLse <- mod$item$se
    print( round(dfr,3)  )
    return( list( "item"=dfr, "jack.itemdiff"=j.itemdiff ))
    }
##########################################################################################
