## File Name: sirt_rmvnorm.R
## File Version: 0.09
## File Last Change: 2018-12-30


sirt_rmvnorm <- function (n, mean=NULL, sigma, ...)
{
    if (is.null(mean)){
        mean <- rep(0,ncol(sigma) )
    }
    CDM::CDM_rmvnorm( n=n, mean=mean, sigma=sigma )
}
