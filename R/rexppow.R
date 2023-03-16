## File Name: rexppow.R
## File Version: 0.22
## File Last Change: 2020-04-27


rexppow <- function (n, mu=0, sigmap=1, pow=2, xbound=100, xdiff=.01)
{
    use_rcpp <- TRUE
    xmin <- -xbound
    xmax <- xbound
    x <- sigmap*( mu + seq(xmin, xmax, by=xdiff) )
    nx <- length(x)
    dx <- dexppow(x=x, mu=mu, sigmap=sigmap, pow=pow)
    dx <- dx / sum(dx)
    dfr <- data.frame( x=x, min=x-xdiff/2, max=x+xdiff/2, cum=cumsum(dx))
    rn <- stats::runif(n)
    if (!use_rcpp){
        res <- rep(NA,n)
        for (nn in 1:n){
            ind_nn <- min( which( rn[nn] <=dfr$cum ) )
            dfr_nn <- dfr[ind_nn,]
            ind_nn1 <- ind_nn-1
            if (ind_nn1==0){
                res[nn] <- dfr_nn[1,"x"]
            } else {
                res[nn] <- dfr_nn$min+xdiff*(rn[nn]-dfr[ind_nn1,"cum"])/
                                    (dfr[ind_nn,"cum"]-dfr[ind_nn1,"cum"])
            }
        }
    } else {
        res <- sirt_rcpp_discrete_inverse(x0=dfr$x, y0=dfr$cum, y=rn)$x
    }
    return(res)
}
