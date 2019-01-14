## File Name: L1_polish.R
## File Version: 0.26



L1_polish <- function(x, eps=0.01, maxiter=30, trace.iter=FALSE, type=1)
{
    #- original median polish procedure
    if (type==1){
        res <- stats::medpolish(x=x, eps=eps, maxiter=maxiter, trace.iter=trace.iter,
                    na.rm=TRUE)
        res <- list(overall=res$overall, row=res$row, col=res$col,
                        residuals=res$residuals)
    }
    #- alternating  median computation
    if (type==2){
        nr <- nrow(x)
        nc <- ncol(x)
        xr <- rep(0,nr)
        xc <- rep(0,nc)
        iterate <- TRUE
        iter <- 1
        while(iterate){
            xc0 <- xc
            xr0 <- xr
            #* update study effects
            delta <- x - xr
            xc <- apply(delta, 2L, stats::median, na.rm=TRUE)
            #* update item parameters
            delta <- x - sirt_matrix2(xc, nrow=nr)
            xr <- apply(delta, 1L, stats::median, na.rm=TRUE)
            change <- max(abs(xc0-xc), abs(xr0-xr), na.rm=TRUE)
            iter <- iter + 1
            if (iter > maxiter){ iterate <- FALSE }
            if (change < eps ){ iterate <- FALSE }
        }
        resid <- x - sirt_matrix2(xc, nrow=nr) - xr
        res <- list(iter=iter, overall=0, row=xr, col=xc, residuals=resid)
    }
    return(res)
}
