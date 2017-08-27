## File Name: mcmc_vcov.R
## File Version: 0.04
## File Last Change: 2017-01-18 11:02:49


###########################################
# variance covariance matrix
mcmc_vcov <- function( mcmcobj , exclude = "deviance" ){	
	mcmcobj <- mcmcobj[ , ! ( colnames(mcmcobj) %in% exclude ) ]
	res <- stats::var(mcmcobj)
	colnames(mcmcobj) -> colnames(res) -> rownames(res)
	return(res)
		}
