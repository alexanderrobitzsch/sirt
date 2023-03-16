## File Name: sirt_moving_average.R
## File Version: 0.14


#----- moving window average for time series
sirt_moving_average <- function(x, B, fill=TRUE)
{
    x1 <- cumsum(x)
    N <- length(x)
    y <- rep(NA,N)
    i <- seq(B+1, N-B)
    xdiff <- x1[ -seq(1,B) ] - x1[ -seq(N-B+1,N) ]
    xdiff <- xdiff[ - seq(1,B) ]
    y[i]  <- ( x1[i] + xdiff - c(0,x1[ -seq(N-2*B,N) ]) ) / (2*B+1)

    # fill NAs at beginning and end of time series
    if(fill){
        j <- seq(0,B-1)
        ybeg <- sapply(j, function(z) sum( x[ seq(1,(2*z+1)) ]) / (2*z+1) )
        yend <- sapply(rev(j), function(z) sum( x[ seq(N-2*z,N) ] ) / (2*z+1) )
        y[j+1] <- ybeg
        y[ rev(N-j) ] <- yend
    }
    return(y)
}


.movingAverage <- sirt_moving_average
