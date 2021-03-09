## File Name: rasch.pairwise.itemcluster.R
## File Version: 0.32


#***** Pairwise estimation with itemclusters
rasch.pairwise.itemcluster <- function( dat, itemcluster=NULL,
            b.fixed=NULL, weights=NULL, conv=.00001, maxiter=3000,
            progress=TRUE, b.init=NULL, zerosum=FALSE)
{
    s1 <- Sys.time()
    CALL <- match.call()
    if ( is.null(b.init) ){
        b.init <- - stats::qlogis( colMeans( dat, na.rm=TRUE ) )
    }
    if (is.null(weights)){
        weights <- rep(1,nrow(dat))
    }
    weights0 <- weights
    I <- ncol(dat)
    dat <- as.matrix(dat)
    dat00 <- dat
    dat[ is.na(dat) ] <- 9
    b <- b.init
    if ( ! is.null(b.fixed) ){
        b[ b.fixed[,1] ] <- b.fixed[,2]
        b.fixed <- cbind( b.fixed, exp( b.fixed[,2] ) )
        zerosum <- FALSE
    }
    # create count tables
    weights <- weights / sum(weights) * nrow(dat)
    sw <- sqrt(weights)
    dat0 <- sw*(dat==0)
    dat1 <- sw*(dat==1)
    Aij <- crossprod(dat0, dat1)
    Aji <- t(Aij)
    # set some entries to zero for itemclusters
    clusters <- unique( itemcluster[ itemcluster !=0 ] )
    CC <- length(clusters)
    for (cc in clusters){
        icc <- which( itemcluster==cc )
        Aji[icc,icc] <- Aij[icc,icc] <- 0
    }
    nij <- Aij + Aji
    Aij_rowsums <- rowSums(Aij)
    eps0 <- eps <- exp(b)
    max.change <- 10
    iter <- 1

    #**** start algorithm
    while( max.change > conv ){
        b0 <- b
        eps0 <- eps
        m1 <- sirt_matrix2( eps0, nrow=I) + matrix( eps0, nrow=I, ncol=I )
        g1 <- rowSums(nij/m1)
        eps <- Aij_rowsums/g1
        b <- log(eps)
        # include item parameter constraints
        if ( ! is.null(b.fixed) ){
            eps[ b.fixed[,1] ] <- b.fixed[,3]
        }
        if (zerosum){
            b1 <- -b
            b2 <- b1-mean(b1)
            eps <- exp(-b2)
            b <- log(eps)
        }
        max.change <- max(abs( b - b0 ))
        if (progress){
            cat( "PL Iter.", iter, ": max. parm. change=", round( max.change, 6 ), "\n")
            utils::flush.console()
        }
        iter <- iter + 1
    }  #** end algorithm

    #** post-processing
    item <- data.frame( N=colSums(1-is.na(dat00)),
                    p=colMeans(dat00, na.rm=TRUE), b=log(eps) )
    if ( is.null(itemcluster) ){
        itemcluster <- rep(0,I)
    }
    item$itemcluster <- itemcluster

    #-- output
    s2 <- Sys.time()
    res <- list( b=b, eps=eps, iter=iter, conv=conv, dat=dat00,item=item,
            fct='rasch.pairwise.itemcluster', itemcluster=itemcluster,
            s1=s1, s2=s2, CALL=CALL )
    class(res) <- 'rasch.pairwise'
    return(res)
}
