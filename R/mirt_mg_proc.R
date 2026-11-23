## File Name: mirt_mg_proc.R
## File Version: 0.116


mirt_mg_proc <- function(modlist)
{
    requireNamespace('mirt')
    requireNamespace('MASS')
    requireNamespace('Matrix')

    G <- length(modlist)
    item <- NULL
    dfr <- NULL
    coef1 <- NULL
    itempars <- NULL
    N <- rep(NA,G)
    logLik <- 0

    for (gg in 1L:G){
        mod_gg <- modlist[[gg]]
        d1 <- mod_gg@Data$data
        N[gg] <- nrow(d1)
        item1 <- colnames(d1)
        I1 <- length(item1)
        var_gg <- mirt::extract.mirt(mod_gg, 'vcov')
        coef_gg <- mirt::extract.mirt(mod_gg, 'parvec')
        logLik <- logLik + mirt::extract.mirt(mod_gg, 'logLik')
        ind1 <- seq(2,2*I1,2)
        ind2 <- seq(1,2*I1,2)
        ind <- c(ind1,ind2)
        coef_gg <- coef_gg[ ind ]
        names(coef_gg) <- c( paste0('G', gg, '_d.',item1), paste0('G', gg, '_a.',item1) )
        c1 <- coef_gg
        coef1 <- c( coef1, coef_gg)
        varmat_gg <- var_gg[ ind, ind ]
        rownames(varmat_gg) <- colnames(varmat_gg) <- names(coef_gg)
        if (gg==1){
            varmat <- varmat_gg
        } else {
            varmat <- Matrix::bdiag(varmat, varmat_gg)
        }

        dfr1 <- data.frame(group=gg, item=item1,
                        par=names(coef_gg), partype=rep(c('d','a'), each=I1),
                        est=coef_gg )
        dfr <- rbind(dfr, dfr1)
        i1 <- data.frame(study=gg, item=item1[1L:I1],
                            a=c1[I1+1L:I1], b=-c1[1L:I1]/c1[I1+1L:I1] )
        itempars <- rbind(itempars, i1)

    }
    rownames(varmat) <- colnames(varmat) <- names(coef1)
    varmat <- as.matrix(varmat)
    infomat <- as.matrix( MASS::ginv(varmat) )
    rownames(infomat) <- colnames(infomat) <- names(coef1)

    # define item identifiers
    dfr$itemid <- match( dfr$item, unique(sort(dfr$item)) )
    I <- max(dfr$itemid)
    NP <- nrow(dfr)/2

    t1 <- table(dfr$itemid)
    i1 <- names(t1)[ t1==2 ]
    dfr$unique <- 0
    dfr[ dfr$itemid %in% i1, 'unique'] <- 1

    ind_a <- which(dfr$partype=='a')
    ind_d <- which(dfr$partype=='d')
    dfr[ ind_a, 'parid'] <- 1L:NP
    dfr[ ind_d, 'parid'] <- 1L:NP

    items <- stats::aggregate(dfr$item, list(dfr$itemid),
                        FUN=function(ll){ ll[1] } )[,2]

    #- include precision weights in itempars
    y <- coef1[ind_a]
    A <- diag(1/y)
    V <- varmat[ind_a, ind_a]
    V <- A %*% V %*% t(A)
    w <- extract_diag(X=V)
    dfr[ind_a, 'weight'] <- itempars$weight_a <- 1/w
    w <- extract_diag(X=varmat, indices=ind_d)
    dfr[ind_d, 'weight'] <- itempars$weight_d <- 1/w

    #-- output
    res <- list( coef=coef1, infomat=infomat, varmat=varmat, parinfo=dfr,
                    itempars=itempars, G=G, I=I, items=items, N=N, logLik=logLik)
    return(res)
}

