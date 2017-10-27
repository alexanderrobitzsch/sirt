## File Name: rm_sdt_pem_inits.R
## File Version: 0.01

rm_sdt_pem_inits <- function( tau.item, a.item, c.rater, d.rater, est.a.item, a.item.fixed, skillspace, PEM )
{
	parmlist <- list( tau.item=tau.item, a.item=a.item, c.rater=c.rater, d.rater=d.rater, mu=0, sigma=1)
	pem_parameter_index <- sirt_pem_create_parameter_index( parmlist=parmlist )
	pem_parameter_sequence <- list()
	center_log_a <- ( ! is.null( a.item.fixed ) ) & ( est.a.item )
	pem_pars <- c("c.rater","d.rater","a.item","tau.item","mu","sigma")	
	if ( skillspace == "discrete" ){
		PEM <- FALSE
	}
	#--- output
	res <- list(PEM=PEM, pem_pars=pem_pars, center_log_a=center_log_a, pem_parameter_index=pem_parameter_index,
				pem_parameter_sequence=pem_parameter_sequence )
	return(res)	
}
