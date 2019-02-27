## File Name: wle.rasch.R
## File Version: 1.294


#-- WLE ability estimation
wle.rasch <- function( dat, dat.resp=NULL, b, itemweights=1+0*b,
            theta=rep(0, nrow(dat)), conv=.001, maxit=200,
            wle.adj=0, progress=FALSE)
{
    theta.change <- 1
    iter <- 0
    N <- nrow(dat)
    if ( is.null(dat.resp) ){
        dat.resp <- 1 - is.na( dat)
        # multiply response matrix with itemweights
        dat.resp <- dat.resp * sirt_matrix2(itemweights, nrow=N)
    }
    dat[ is.na( dat) ] <- 9
    sufftheta <- rowSums( dat.resp *  dat  )
    if (wle.adj>0){
        nr <- rowSums( dat.resp )
        i1 <- which( nr==sufftheta )
        if ( length(i1) > 0 ){
            sufftheta[ i1 ] <- nr[ i1 ] - wle.adj
        }
        i1 <- which( sufftheta==0 )
        if ( length(i1) > 0 ){
            sufftheta[ i1 ] <- wle.adj
        }
    }
    if ( progress){
        cat("\n  WLE estimation  |" )
    }
    old_increment <- rep( 5, length(theta))
    eps <- 1E-10
    while( ( max(abs( theta.change )) > conv ) & ( iter < maxit ) ){
        # calculate P and Q
        p.ia <- stats::plogis( theta - sirt_matrix2(b, nrow=N) )
        q.ia <- 1 - p.ia
        # Log Likelihood (for every subject)
        M0 <- dat.resp * p.ia
        l1 <- sufftheta - rowSums(M0)
        # I and J terms
        M1 <- M0*q.ia
        I.wle <- rowSums(M1)
        J.wle <- rowSums(M1*(q.ia - p.ia))
        I1.wle <- J.wle
        J1.wle <- rowSums( M1*( 6 * p.ia^2 - 6 * p.ia + 1 ) )
        # calculate objective function
        f.obj <- l1 + J.wle / ( 2 * I.wle + eps)
        # derivative of the objective function
        f1.obj <- - I.wle + ( J1.wle * I.wle - J.wle  * I1.wle )/ ( 2 * I.wle^2 + eps)
        # theta change
        theta.change <- - f.obj / ( f1.obj + eps )
        # define damped increment!!
        max_increment <- abs(old_increment)
        theta.change <- sirt_trim_increment(increment=theta.change, max_increment=max_increment)
        old_increment <- abs(theta.change)
        theta <- theta + theta.change
        iter <- iter + 1
        if ( any( is.nan( theta.change ) ) ){
            stop( "Numerical problems occur during WLE estimation procedure.")
        }
        if ( progress){
            cat("-")
        }
    }
    res <- list( theta=theta, se.theta=1/sqrt(abs(f1.obj)),
                        dat.resp=dat.resp, p.ia=p.ia )
    #*** compute WLE reliability
    v1 <- stats::var(res$theta)
    v2 <- mean(res$se.theta^2 )
    wle.rel <- ( v1 - v2 ) / v1
    cat("\nWLE Reliability=", round(wle.rel,3), "\n")
    res$wle.rel <- wle.rel
    return(res)
}
