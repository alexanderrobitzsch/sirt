## File Name: mcmc_extract_samples_first_chain.R
## File Version: 0.04
## File Last Change: 2018-12-30

mcmc_extract_samples_first_chain <- function(mcmcobj)
{
    if ( ! is.matrix(mcmcobj) ){
        dat.bugs <- mcmcobj[[1]]
    } else {
        dat.bugs <- mcmcobj
    }
    return(dat.bugs)
}
