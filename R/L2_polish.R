## File Name: L2_polish.R
## File Version: 0.07
## File Last Change: 2019-01-14

L2_polish <- function(x, eps=1e-4, maxiter=20)
{
    x0 <- x
    nr <- nrow(x)
    nc <- ncol(x)
    wgt <- 1*(! is.na(x))
    iterate <- TRUE
    xr <- 0
    xc <- 0
    iter <- 1
    while(iterate){
        xr_old <- xr
        xc_old <- xc
        x01 <- x0 - xr
        xc <- weighted_colMeans( mat=x01, wgt=wgt )
        xc <- xc - mean(xc)
        x01 <- x0 - sirt_matrix2(xc, nrow=nr)
        xr <- weighted_rowMeans( mat=x01, wgt=wgt )
        change <- max( abs(xc_old-xc), abs(xr_old-xr), na.rm=TRUE)
        if (iter > maxiter){ iterate <- FALSE }
        if (change < eps ){ iterate <- FALSE }
    }
    #--- output
    res <- list(row=xr, col=xc, iter=iter, wgt=wgt)
    return(res)
}
