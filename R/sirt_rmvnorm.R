## File Name: sirt_rmvnorm.R
## File Version: 0.04
## File Last Change: 2017-09-19 20:46:42


sirt_rmvnorm <- function (n, mean=NULL, sigma, ...) 
{
	if (is.null(mean)){
		mean <- rep(0,ncol(sigma) )
	}	
	CDM::CDM_rmvnorm( n=n, mean=mean, sigma = sigma )
}
