## File Name: ccov_np_regression.R
## File Version: 0.07

ccov_np_regression <- function(x, y, xgrid, bwscale=1.1, smooth=TRUE, score=NULL)
{
    N <- length(x)
    if (smooth){
        y <- stats::ksmooth( x=x, y=y, bandwidth=bwscale*N^(-1/5),
                            x.points=xgrid, kernel="normal")$y
    } else {
        a1 <- stats::aggregate(y, list(score), mean, na.rm=TRUE)
        y <- a1[,2]
    }
    return(y)
}
