## File Name: wle.rasch.R
## File Version: 1.20

#----------------------------------------------------------
# Function for WLE ability estimation
wle.rasch <- function( dat, dat.resp=NULL, b, # a=1 + 0*b, c=0*b,
                            itemweights=1+0*b,
                            theta=rep(0, nrow(dat)), conv=.001, maxit=200,
                            wle.adj=0, progress=FALSE)
{
    #-------------------------------------------------------#
    # INPUT:                                                #
    # dat   ... data frame with item response pattern       #
    # dat.resp ... response pattern (non-missing items)     #
    # theta ... initial values for ability estimation       #
    # itemweights itemweights for likelihood estimation     #
    # b     ... item difficulty parameters                  #
    # a     ... item discriminations                        #
    # c     ... guessing parameter                          #
    # conv  ... convergence criterion                       #
    #-------------------------------------------------------#
    theta.change <- 1
    iter <- 0
#    if ( sum(theta)==0 ){
#        theta <- qnorm( (rowMeans( dat, na.rm=T ) + .01 ) / 1.02 )
#                }
    if ( is.null(dat.resp) ){
            dat.resp <- 1 - is.na( dat)
            # multiply response matrix with itemweights
            dat.resp <- dat.resp * outer( rep(1,nrow(dat)), itemweights )
                            }
            dat[ is.na( dat) ] <- 9
    sufftheta <- rowSums( dat.resp *  dat  )
    if (wle.adj>0){
        nr <- rowSums( dat.resp )
        i1 <- which( nr==sufftheta )
        if ( length(i1) > 0 ){    sufftheta[ i1 ] <- nr[ i1 ] - wle.adj }
        i1 <- which( sufftheta==0 )
        if ( length(i1) > 0 ){    sufftheta[ i1 ] <- wle.adj }
                    }
    if ( progress){
        cat("\n  WLE estimation  |" )
    }
    old_increment <- rep( 5, length(theta))
    while( ( max(abs( theta.change )) > conv ) & ( iter < maxit ) ){
        # calculate P and Q
        p.ia <- stats::plogis( outer( theta, b, "-" ) )
        q.ia <- 1 - p.ia
        # Log Likelihood (for every subject)
        l1 <- sufftheta - rowSums( dat.resp * p.ia )
        # I and J terms
        M1 <- dat.resp * p.ia * q.ia
        I.wle <- rowSums( M1 )
        I.wle <- rowSums( dat.resp * p.ia * q.ia )
        J.wle <- rowSums( M1 * (q.ia - p.ia ) )
        I1.wle <- J.wle
        J1.wle <- rowSums( M1* ( 6 * p.ia^2 - 6 * p.ia + 1 ) )
        # calculate objective function
        f.obj <- l1 + J.wle / ( 2 * I.wle )
        # derivative of the objective function
        f1.obj <- - I.wle + ( J1.wle * I.wle - J.wle  * I1.wle )/ ( 2 * I.wle^2 )
        # theta change
        increment <- theta.change <- - f.obj / f1.obj
        # define damped increment!!
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 1E-10 ) )
        theta.change <- ifelse( abs( increment) > abs(old_increment), increment/(2*ci), increment )
        old_increment <- abs(theta.change)
        theta <- theta + theta.change
        iter <- iter + 1
        if ( any( is.nan( theta.change ) ) ){
                stop( "Numerical problems occur during WLE estimation procedure.")
                        }
        if ( progress){  cat("-") }
    }
    res <- list( "theta"=theta, "se.theta"=1 / sqrt( abs( f1.obj ) ),
                "dat.resp"=dat.resp,    "p.ia"=p.ia )
    #*** compute WLE reliability
    v1 <- stats::var(res$theta)
    v2 <- mean(res$se.theta^2 )
    wle.rel <- ( v1 - v2 ) / v1
    cat("WLE Reliability=", round(wle.rel,3), "\n")
    res$wle.rel <- wle.rel
    return(res)
 }
#--------------------------------------------------------------------------------------------------------#



