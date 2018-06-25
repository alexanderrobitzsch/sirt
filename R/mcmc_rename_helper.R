## File Name: mcmc_rename_helper.R
## File Version: 0.02

mcmc_rename_helper <- function( string, rep_string=3, pre=3, suff=3)
{
    string <- paste0(rep( string, rep_string ), collapse="")
    pre <- paste0(rep( "_", pre ), collapse="")
    suff <- paste0(rep( "_", suff ), collapse="")
    trans <- paste0( pre, string, suff )
    return(trans)
}
