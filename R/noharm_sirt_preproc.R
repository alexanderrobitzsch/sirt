## File Name: noharm_sirt_preproc.R
## File Version: 0.409


#**** data preprocessing noharm.sirt
noharm_sirt_preproc <- function( dat, pm, N, weights, Fpatt, Fval,
    Ppatt, Pval, Psipatt, Psival, wgtm, dimensions, pos.loading,
    pos.variance, pos.residcorr, input_pm, lower, upper )
{

    res <- NULL
    res$dat0 <- dat

    if ( ! input_pm ){
        N <- nrow(dat)
        I <- ncol(dat)
        items <- colnames(dat)
    } else {
        I <- ncol(pm)
        items <- colnames(pm)
        pm <- as.matrix(pm)
    }

    if (length(lower)==1){
        lower <- rep(lower, I)
    }
    if (length(upper)==1){
        upper <- rep(upper, I)
    }

    #---- data processing
    res$N <- N
    res$I <- I

    if (! input_pm){
        if (is.null(weights)){
            weights <- rep(1, N)
        }
        res$weights <- weights
        dat.resp <- 1-is.na(dat)
        dat0 <- dat
        dat[dat.resp==0] <- 0
        # calculate (weighted) product moment correlation
        ss <- as.matrix( crossprod( weights*dat.resp ) )
        eps <- 1e-6
        pm <- as.matrix( crossprod( as.matrix(dat*weights*dat.resp) ) / ( ss+eps ) )
    } else {
        dat.resp <- NULL
        ss <- 1
    }
    res$dat <- dat
    res$dat.resp <- dat.resp
    res$pm <- pm
    res$pm0 <- pm
    res$ss <- ss
    res$lower <- lower
    res$upper <- upper

    # CFA or EFA?
    if ( is.null(dimensions) ){
        model.type <- 'CFA'
        modtype <- 3    # 3 - multidimensional CFA
    } else {
        model.type <- 'EFA'
        modtype <- 2    # 2 - multidimensional EFA
        D <- dimensions
        Pval <- diag(D)
        Ppatt <- 0*diag(D)
        Fpatt <- matrix(1,nrow=I,ncol=D)
        if (D>1){
            for (dd in 2L:D){
                Fpatt[dd,1L:(dd-1)] <- 0
            }
        }
        Fval <- .5*(Fpatt>0)
        if ( D==1 ){ # 1 dimension
            model.type <- 'CFA'
            modtype <- 3
        }
    }

    res$model.type <- model.type
    res$modtype  <- modtype
    # initial values if they are not provided
    if ( is.null(Psival) ){ Psival <- 0*diag(res$I) }
    if ( is.null(Psipatt) ){ Psipatt <- 0*diag(res$I) }
    if ( is.null(Fval) ){ Fval <- .5*(Fpatt>0) }
    if ( is.null(Pval) ){ Pval <- diag( ncol(Ppatt) ) }
    # weight matrix
    wgtm.default <- FALSE
    if ( is.null(wgtm) ){
        wgtm <- matrix(1, nrow=I, ncol=I)
        wgtm.default <- TRUE
    }
    diag(wgtm) <- 0
    wgtm <- ( wgtm + t(wgtm) ) / 2
    wgtm <- wgtm * ( ss > 0 )
    res$wgtm <- wgtm
    res$sumwgtm <- ( sum( wgtm > 0 ) - sum( diag(wgtm) > 0 ) ) / 2

    #*** column names
    D <- ncol(Ppatt)
    cn <- paste0('F',1L:D)
    if (is.null(colnames(Fpatt) ) ){
        colnames(Fpatt) <- cn
    }
    if (is.null(colnames(Fval) ) ){
        colnames(Fval) <- colnames(Fpatt)
    }

    if (is.null(colnames(Pval))){
        colnames(Pval) <- colnames(Ppatt)
    }
    if (is.null(colnames(Pval))){
        colnames(Pval) <- colnames(Fval)
    }
    rownames(Pval) <- colnames(Pval)

    #--- create parameter table

    # F
    parm_table <- noharm_sirt_preproc_parameter_table_matrix(pattmat=Fpatt, valmat=Fval,
            patt_id=1, patt_label='F', minval=0, symm=FALSE)
    # P
    parm1 <- noharm_sirt_preproc_parameter_table_matrix(pattmat=Ppatt, valmat=Pval,
            patt_id=2, patt_label='P', minval=sirt_max(parm_table$index), symm=TRUE)
    parm_table <- rbind(parm_table, parm1)

    # Psi
    parm1 <- noharm_sirt_preproc_parameter_table_matrix(pattmat=Psipatt, valmat=Psival,
            patt_id=3, patt_label='Psi', minval=sirt_max(parm_table$index), symm=TRUE)
    parm_table <- rbind(parm_table, parm1)
    parm_table <- parm_table[ parm_table$nonnull_par==1, ]
    rownames(parm_table) <- NULL
    npar <- max(parm_table$index, na.rm=TRUE)

    # indices
    ip <- parm_table$index
    ip[duplicated(ip)] <- NA
    extract_index <- match(1L:npar, ip)
    extract_index <- intersect( which( ! is.na( parm_table$index ) ),
                            which( ! duplicated( parm_table$index ) ) )
    parm_table$est <- parm_table$starts

    parm_table$lower <- -Inf

    ind <- NULL
    if (pos.variance){
        ind1 <- which((parm_table$mat=='P') & (parm_table$row==parm_table$col ))
        ind <- union(ind, ind1)
    }
    if (pos.loading){
        ind1 <- which((parm_table$mat=='F') )
        ind <- union(ind, ind1)
    }
    if (pos.residcorr){
        ind1 <- which((parm_table$mat=='Psi') )
        ind <- union(ind, ind1)
    }
    parm_table[ind, 'lower'] <- 0
    parm_table[ is.na(parm_table$index), 'lower'] <- NA
    non_fixed <- ! parm_table$fixed
    include_index <- parm_table$index[ non_fixed ]
    parm_table$nonnull_par <- NULL
    attr(parm_table, 'extract_index') <- extract_index
    attr(parm_table, 'non_fixed') <- non_fixed
    attr(parm_table, 'include_index') <- include_index
    attr(parm_table, 'npar') <- npar
    attr(parm_table, 'NH') <- sum(non_fixed)
    attr(parm_table, 'est_par_index') <- which(parm_table$est_par==1)
    attr(parm_table, 'parm_table_free_index') <- which(parm_table$fixed==0)

    parm_index <- list()
    for (mat_label in c('F', 'P', 'Psi') ){
        ind_mat <- which(parm_table$mat==mat_label)
        parm_index[[ mat_label ]] <- list()
        parm_index[[ mat_label ]][[ 'row_parm_table' ]] <- ind_mat
        parm_index[[ mat_label ]][[ 'entries' ]] <- parm_table[ ind_mat, c('row','col')]
        parm_index[[ mat_label ]][[ 'len' ]] <- length(ind_mat)
        nrow <- I
        ncol <- D
        if (mat_label=='P'){ nrow <- D}
        if (mat_label=='Psi'){ ncol <- I}
        parm_index[[ mat_label ]][['nrow']] <- nrow
        parm_index[[ mat_label ]][['ncol']] <- ncol
    }

    res$parm_table <- parm_table
    res$npar <- npar
    res$parm_index <- parm_index

    #*****
    # matrix conversion
    #****
    res$Fpatt <- as.matrix(Fpatt)
    res$Fval <- as.matrix(Fval)
    res$Ppatt <- as.matrix(Ppatt)
    res$Pval <- as.matrix(Pval)
    res$Psipatt <- as.matrix(Psipatt)
    res$Psival <- as.matrix(Psival)
    res$Psipatt <- ( Psipatt + t(Psipatt) )/2
    res$Ppatt <- ( Ppatt + t(Ppatt) )/2
    res$D <- ncol(Ppatt)
    res$Fpatt <- 1*(res$Fpatt>0)
    res$Ppatt <- 1*(res$Ppatt>0)
    res$Psipatt <- 1*(res$Psipatt>0)
    res$wgtm.default <- wgtm.default
    res$F_dimnames <- colnames(Fpatt)
    res$items <- items
    return(res)
}


