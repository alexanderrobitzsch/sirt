## File Name: mcmc_extract_samples_first_chain.R
## File Version: 0.02

mcmc_extract_samples_first_chain <- function(mcmcobj)
{
    if ( ! is.matrix(mcmcobj) ){
        dat.bugs <- mcmcobj[[1]]
    } else {
        dat.bugs <- mcmcobj
    }
    return(dat.bugs)
}
