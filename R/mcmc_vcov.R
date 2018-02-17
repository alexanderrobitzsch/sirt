## File Name: mcmc_vcov.R
## File Version: 0.06


###########################################
# variance covariance matrix
mcmc_vcov <- function( mcmcobj , exclude = "deviance" )
{	
	mcmcobj <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
	mcmcobj <- mcmcobj[ , ! ( colnames(mcmcobj) %in% exclude ) ]
	res <- stats::var(mcmcobj)
	colnames(mcmcobj) -> colnames(res) -> rownames(res)
	return(res)
}
