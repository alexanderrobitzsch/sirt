## File Name: mcmc_coef.R
## File Version: 0.09

###########################################
# coefficients from one MCMC chain
mcmc_coef <- function( mcmcobj, exclude="deviance" )
{
    mcmcobj <- mcmc_extract_samples_first_chain(mcmcobj=mcmcobj)
    mcmcobj <- mcmcobj[, ! ( colnames(mcmcobj) %in% exclude ) ]
    res <- colMeans(mcmcobj)
    colnames(mcmcobj) -> names(res)
    return(res)
}
