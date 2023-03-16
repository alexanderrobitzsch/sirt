## File Name: rasch.pairwise.R
## File Version: 0.446


#------ Rasch estimation (Approximate Method)
# also called MINCHI method
# Handbook of Statistics Vol. 26
# Chapter of G. Fischer: p. 544 ff.
rasch.pairwise <- function( dat, weights=NULL, conv=0.0001, maxiter=3000,
        progress=TRUE, b.init=NULL, zerosum=FALSE, power=1,
        direct_optim=TRUE)
{
    s1 <- Sys.time()
    CALL <- match.call()
    version <- 1
    # should items be excluded?
    item.elim <- which( colMeans( dat, na.rm=TRUE ) %in% c(0,1))
    if (length(item.elim)>0){
        dat <- dat[, - item.elim ]
    }
    I <- ncol(dat)
    dat <- as.matrix(dat)
    dat00 <- dat
    if (direct_optim){
        p <- colMeans(dat, na.rm=TRUE)
        p <- p[ order(abs(p-.5)) ]
        ind0 <- match(names(p),colnames(dat))
        dat <- dat[,ind0]
    }

    # data preparation
    dat.resp <- 1 - is.na(dat)
    dat.9 <- dat
    dat.9[ is.na(dat) ] <- 9
    # calculate n_{ij}
    if (is.null(weights)){
        weights <- rep(1,nrow(dat))
    }
    sw <- sqrt(weights)
    n.ij <- crossprod( dat.9 * dat.resp * sw, ( 1 - dat.9 ) * dat.resp *sw )
    # which item pairs occur in estimation procedure?
    delta.ij <- 1 * ( n.ij + t(n.ij) > 0 )

    # initial values for beta
    if( is.null(b.init) ){
        beta <- - stats::qlogis( colMeans( dat, na.rm=TRUE ) )
    } else {
        beta <- b.init
    }
    # calculate y_{ij} values
    if (power==1){
        y.ij <- n.ij / ( n.ij + t(n.ij) )
    }
    if (power>1){
        m.ij <- n.ij
        for (hh in 2L:power){
            m.ij <- m.ij %*% n.ij
        }
        y.ij <- m.ij / ( n.ij + t(n.ij) )  # Fischer (2007)
    }
    if (version==2){
        y.ij <- n.ij^2 / ( n.ij + t(n.ij) )
    }
    y.ij[ delta.ij==0 ] <- 0
    y.ji <- t(y.ij)
    eps <- exp(-beta)
    b <- -log(eps)

    if (direct_optim){
        #- direct optimization
        n.ji <- t(n.ij)
        res <- rasch_pairwise_optimize(n.ij=n.ij, n.ji=n.ji, beta=beta,
                    zerosum=zerosum, maxiter=maxiter, estimator="MINCHI")
    } else {
        #- iterations
        fix_first <- version==2
        res <- rasch_pairwise_iterations( eps=eps, y.ij=y.ij, delta.ij=delta.ij,
                    conv=conv, maxiter=maxiter, zerosum=zerosum, progress=progress,
                    fix_first=fix_first)
    }
    eps <- res$eps
    b <- res$b
    iter <- res$iter

    #* post processing
    item <- data.frame( item=colnames(dat), N=colSums(1-is.na(dat)),
                            p=colMeans(dat, na.rm=TRUE), b=b, eps=eps,
                            itemcluster=rep(0,I) )
    item <- item[ match(colnames(dat00), item$item), ]

    #-- output
    s2 <- Sys.time()
    res <- list( b=item$b, eps=item$eps, iter=iter, conv=conv, dat=dat, I=I,
                freq.ij=n.ij, item=item, power=power, direct_optim=direct_optim,
                fct='rasch.pairwise', s1=s1, s2=s2, CALL=CALL )
    class(res) <- 'rasch.pairwise'
    return(res)
}

