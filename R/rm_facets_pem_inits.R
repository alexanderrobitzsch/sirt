## File Name: rm_facets_pem_inits.R
## File Version: 0.06

rm_facets_pem_inits <- function( tau.item, a.item, a.rater, b.rater, skillspace, PEM,
        a.item.fixed, est.a.item)
{
    parmlist <- list( tau.item=tau.item, a.item=a.item, b.rater=b.rater, a.rater=a.rater,
                        mu=0, sigma=1)
    pem_parameter_index <- sirt_pem_create_parameter_index( parmlist=parmlist )
    pem_parameter_sequence <- list()
    center_log_a <- ( ! is.null( a.item.fixed ) ) & ( est.a.item )
    pem_pars <- c("b.rater","a.rater","a.item","tau.item","mu","sigma")
    if ( skillspace=="discrete" ){
        PEM <- FALSE
    }
    #--- output
    res <- list(PEM=PEM, pem_pars=pem_pars, center_log_a=center_log_a,
                    pem_parameter_index=pem_parameter_index,
                    pem_parameter_sequence=pem_parameter_sequence )
    return(res)
}
