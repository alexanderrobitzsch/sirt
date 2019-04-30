## File Name: polychoric2.R
## File Version: 0.328


#---- estimating polychoric correlation using the Olsson method
# of maximum likelihood estimation
polychoric2 <- function( dat, maxiter=100, cor.smooth=TRUE, use_pbv=TRUE,
    conv=1e-10, rho_init=NULL, weights=NULL )
{
    dat1 <- as.matrix(dat)
    NV <- ncol(dat1)
    N <- nrow(dat1)
    # initial values for rho
    if (is.null(rho_init)){
        rho_init <- matrix(0, nrow=NV, ncol=NV)
    }
    if (is.null(weights)){
        weights <- rep(1,N)
    }
    min_val <- apply(dat1, 2, min, na.rm=TRUE)
    if (any(min_val>0)){
        stop("Minimum value must always zero.\n")
    }

    # maximum number of categories
    maxK <- max(dat1, na.rm=TRUE )
    # compute polychoric correlation
    res0 <- sirt_rcpp_polychoric2( dat=dat1, maxK=maxK, maxiter=maxiter,
                    use_pbv=use_pbv, conv=conv, rho_init=rho_init,
                    weights=weights)
    iter <- res0$iter
    rho <- res0$rho
    Nobs <- res0$Nobs
    maxcat <- res0$maxcat
    thresh <- res0$thresh

    #--- output cleaning
    # thresholds
    tau <- thresh[,c(2:(maxK+1))]
    tau[ tau==99 ] <- Inf
    rownames(rho) <- colnames(rho) <- colnames(dat1)
    if ( maxK > 1 ){
        rownames(tau) <- rownames(rho)
        colnames(tau) <- paste0("Cat", 1:maxK)
    }
    if ( maxK==1 ){
        names(tau) <- rownames(rho)
    }
    # handling missing entries in rho
    rho[ Nobs==0 ] <- NA
    diag(rho) <- 1
    if ( sum( is.na(rho) > 0 ) ){
        cor.smooth <- FALSE
    }
    if (cor.smooth){
        rho <- sirt_import_psych_cor.smooth(x=rho)
    }
    # output list
    res <- list(tau=tau, rho=rho, Nobs=Nobs, maxcat=maxcat,
                cor.smooth=cor.smooth, iter=iter)
    return(res)
}

