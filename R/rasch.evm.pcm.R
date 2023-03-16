## File Name: rasch.evm.pcm.R
## File Version: 1.59
## File Last Change: 2018-12-30


# estimation of the partial credit model employing the eigenvector method
rasch.evm.pcm <- function( dat, jackunits=20, weights=NULL,
    pid=NULL, group=NULL, powB=2, adj_eps=.3, progress=TRUE )
{
    CALL <- match.call()
    s1 <- Sys.time()
    dat0 <- dat <- as.matrix(dat)
    N <- nrow(dat)
    if ( is.null(weights) ){
        weights <- rep(1,N)
    }
    dat.resp <- 1-is.na(dat)
    dat[ is.na(dat) ] <- 0
    if ( length(jackunits)==1 ){
        jackunits <- ( 1:N ) %% jackunits
    }
    JJ <- length(unique(jackunits))
    maxK <- sirt_colMaxs(x=dat)
    I <- ncol(dat)
    # row indices of B
    dfr <- NULL
    for (ii in 1:I){
        dfr <- rbind( dfr, cbind( ii-1, 1:maxK[ii] )    )
    }
    row_index <- as.matrix(dfr)
    # column indices of B
    dfr <- NULL
    for (ii in 1:I){
        dfr <- rbind( dfr, cbind( ii-1, 1:maxK[ii] - 1) )
    }
    col_index <- as.matrix(dfr)
    powD <- powB
    if (progress){ progress <- 1 }

    #**************
    # subroutine multiple groups
    nogroup <- TRUE
    G <- 1
    if ( ! is.null(group) ){
        group0 <- group
        group.unique <- sort(unique(group))
        group <- match( group, group.unique)
        G <- length(group.unique)
        nogroup <- FALSE
    }

    #*** descriptives
    if (nogroup){
        group.unique <- NULL
        desc <- data.frame( Nobs=nrow(dat), G=1 )
        desc$N.items <- ncol(dat)
        desc$N.pars <- nrow(dfr)
        desc$sum.weights <- sum(weights)
    }
    if (!nogroup){
        desc <- data.frame( Nobs=as.numeric(table(group)), G=G )
        desc$N.items <- ncol(dat)
        desc$N.pars <- nrow(dfr)
        desc$sum.weights <- stats::aggregate( weights, list(group), sum)[,2]
    }
    m1 <- rowSums(dat.resp)
    desc$M.Nitems <- mean(m1)
    desc$SD.Nitems <- stats::sd(m1)
    desc$min.Nitems <- min(m1)
    desc$max.Nitems <- max(m1)

    #**** run Rcpp subroutine
    # no group
    if ( nogroup ){
        res1 <- sirt_rcpp_evm_comp_poly( dat=dat, dat_resp=dat.resp,
                    weights=weights, JJ=JJ, jackunits=jackunits, powD=powD, progress_=progress,
                    row_index=row_index, col_index=col_index )
    }  else {
        # groups
        res1 <- as.list(1:G)
        for (gg in 1:G){
            if (progress){
                cat("\n------- Group", group.unique[gg], " ------\n")
            }
            ind.gg <- which(group==gg)
            res1[[gg]] <- sirt_rcpp_evm_comp_poly( dat=dat[ind.gg,],
                                dat_resp=dat.resp[ind.gg,], weights=weights[ind.gg], JJ=JJ,
                                jackunits=jackunits, powD=powD, progress_=progress, row_index=row_index,
                                col_index=col_index )
        }
    }
    #************************
    # collect item parameters
    item_index <- row_index[,1]+1
    items <- colnames(dat)[ item_index ]
    ri <- paste0( items, "_Cat", row_index[,2]  )
    # item parameters
    if (nogroup){
        PP <- nrow(res1$PARS_vcov) - 2
    } else {
        PP <- nrow(res1[[1]]$PARS_vcov) - 2
    }
    item <- data.frame(parmlabel=ri, item=item_index, itemlabel=items,
                category=row_index[,2] )
    item$parmindex <- seq(1,PP)
    if ( nogroup){
        b_evm <- res1$b_evm
        item$est <- res1$b_evm
        item$se <- sqrt( diag( res1$PARS_vcov )[1:PP] )
    } else {
        for (gg in 1:G){
            rgg <- res1[[gg]]
            item[,paste0("est.Gr",group.unique[gg])] <- rgg$b_evm
            item[,paste0("se.Gr",group.unique[gg])] <- sqrt( diag( rgg$PARS_vcov )[1:PP] )
        }
    }
    if (nogroup){
        colnames(res1$B) <- rownames(res1$B) <- ri
        # jackknife bias corrected estimation
            # http://statweb.stanford.edu/~susan/courses/s208/node16.html
        item$est_jack <- item$est - (res1$JJadj - 1 ) * ( res1$PARS_means[1:PP] - item$est )
        PARS_vcov <- res1$PARS_vcov[ 1:PP, 1:PP ]
    }
    rownames(item) <- NULL
    #*************************
    # include person parameter estimation
    #   person parameter estimation: no multiple groups
    #   use function mle.pcm.group
    person <- NULL
    if ( nogroup){   # no groups
        row_index[,1] <- row_index[,1] + 1
        b <- matrix( NA, nrow=I, ncol=max( maxK) )
        rownames(b) <- colnames(dat)
        colnames(b) <- paste0("Cat", 1:max(maxK) )
        b[ row_index ] <- item$est
        #****
        # estimation subroutine
        person <- mle.pcm.group(dat, b, a=rep(1, ncol(dat)), pid=pid,
                        adj_eps=adj_eps, conv=1e-04, maxiter=30)$person
        B <- res1$B
        D <- res1$D
        #  PARS_vcov <- PARS_vcov
        JJadj <- res1$JJadj
    }
    if ( ! nogroup){   # groups
        row_index[,1] <- row_index[,1] + 1
        b0 <- as.list(1:G)
        for (gg in 1:G){
            b <- matrix( NA, nrow=I, ncol=max( maxK) )
            rownames(b) <- colnames(dat)
            colnames(b) <- paste0("Cat", 1:max(maxK) )
            b[ row_index ] <- item[,paste0("est.Gr",group.unique[gg])]
            b0[[gg]] <- b
        }
        b <- b0
        b_evm <- B <- D <- PARS_vcov <- JJadj <- as.list(1:G)
        for (gg in 1:G){
            rgg <- res1[[gg]]
            B[[gg]] <- rgg$B
            D[[gg]] <- rgg$D
            PARS_vcov[[gg]] <- rgg$PARS_vcov[ 1:PP, 1:PP ]
            JJadj[[gg]] <- rgg$JJadj
            b_evm[[gg]] <- rgg$b_evm
        }
    }

    #*** assess differential item functioning
    difstats <- NULL
    if ( ! nogroup ){
        difstats <- rasch_evm_pcm_dif( b_evm=b_evm, item=item,
                        PARS_vcov=PARS_vcov, I=I, G=G, group.unique=group.unique, dat=dat,
                        dat.resp=dat.resp )
    }

    # collect coefficients in case one group and multiple groups
    # create results for multiple groups
    if (G>1){
        ni <- nrow(item)
        Npars <- ni*G
        coeff <- rep( NA, Npars )
        pcov <- matrix( 0, nrow=Npars, ncol=Npars )
        for (gg in 1:G){
            ind.gg <- (gg-1)*ni + 1:ni
            coeff[ind.gg] <- item[, paste0("est.Gr", group.unique[gg] ) ]
            names(coeff)[ind.gg] <- paste0( item$parmlabel, "_Gr", group.unique[gg] )
            pcov[ind.gg,ind.gg] <- PARS_vcov[[gg]]
        }
        rownames(pcov) <- colnames(pcov) <- names(coeff)
        PARS_vcov <- pcov
    } else {
        coeff <- as.vector(item$est)
        names(coeff) <- item$parmlabel
        rownames(PARS_vcov) <- colnames(PARS_vcov) <- names(coeff)
    }
    #************************* output *******************************
    s2 <- Sys.time()
    res <- list( item=item, b=b, person=person, B=B, D=D, coef=coeff, vcov=PARS_vcov,
                JJ=JJ, JJadj=JJadj, powB=powB, maxK=maxK, G=G, desc=desc,
                difstats=difstats, b_evm=b_evm, I=I, group.unique=group.unique,
                dat=dat0, CALL=CALL, s1=s1, s2=s2)
    class(res) <- 'rasch.evm.pcm'
    return(res)
}
