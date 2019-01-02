## File Name: rasch.pairwise.R
## File Version: 0.38


#------ Rasch estimation (Approximate Method)
# also called MINCHI method
# Handbook of Statistics Vol. 26
# Chapter of G. Fischer: p. 544 ff.
# pairwise likelihood method
rasch.pairwise <- function( dat, conv=0.0001, maxiter=3000,
        progress=TRUE, b.init=NULL, zerosum=FALSE )
{
    s1 <- Sys.time()
    CALL <- match.call()
    # should items be excluded?
    item.elim <- which( colMeans( dat, na.rm=TRUE ) %in% c(0,1))
    if (length(item.elim)>0){
        dat <- dat[, - item.elim ]
    }

    dat <- as.matrix(dat)
    I <- ncol(dat)
    # data preparation
    dat.resp <- 1 - is.na( dat )
    dat.9 <- dat
    dat.9[ is.na(dat) ] <- 9
    # calculate n_{ij}
    n.ij <- crossprod( dat.9 * dat.resp, ( 1 - dat.9 ) * dat.resp  )
    # which item pairs occur in estimation procedure
    delta.ij <- 1 * ( n.ij + t( n.ij ) > 0 )

    # initial values for beta
    if( is.null( b.init) ){
        beta <- - stats::qlogis( colMeans( dat, na.rm=TRUE ) )
    } else {
        beta <- b.init
    }
    # calculate y_{ij} values
    y.ij <- n.ij / ( n.ij + t( n.ij) )
    y.ij[ delta.ij==0 ] <- 0
    y.ji <- t( y.ij )
    eps <- exp(-beta)
    change <- 1
    iter <- 0
    b <- -log(eps)

    #* start estimation algorithm
    while( change > conv & iter < maxiter ){
        eps0 <- eps
        b0 <- b
        eps <- sqrt( rowSums( y.ij * eps * delta.ij ) / colSums( y.ij / eps ) )
        if (zerosum){
            b1 <- - log(eps)
            b2 <- b1 - mean(b1)
            eps <- exp(-b2)
        }
        b <- -log(eps)
        change <- max( abs( eps0 - eps ) )
        change_b <- max(abs( -b0 + b))
        iter <- iter + 1
        if ( progress ){
            cat( "PL Iter.", iter, ": max. parm. change=", round( change_b, 6 ), "\n")
            utils::flush.console()
        }
    } #* end estimation algorithm

    #* post-processing
    item <- data.frame( N=colSums(1-is.na(dat)), p=colMeans(dat, na.rm=TRUE ),
                        b=b, itemcluster=rep(0,I) )
    #-- output
    s2 <- Sys.time()
    res <- list( b=- log(eps), eps=eps, iter=iter, conv=conv, dat=dat, I=I,
                freq.ij=n.ij, item=item, fct='rasch.pairwise',
                s1=s1, s2=s2, CALL=CALL )
    class(res) <- 'rasch.pairwise'
    return(res)
}

