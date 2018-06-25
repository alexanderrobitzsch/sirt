## File Name: sirt_rmvnorm.R
## File Version: 0.07


sirt_rmvnorm <- function (n, mean=NULL, sigma, ...)
{
    if (is.null(mean)){
        mean <- rep(0,ncol(sigma) )
    }
    CDM::CDM_rmvnorm( n=n, mean=mean, sigma=sigma )
}
