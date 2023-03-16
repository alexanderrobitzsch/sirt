## File Name: mcmc_as_formula.R
## File Version: 0.11
## File Last Change: 2018-12-30

mcmc_as_formula <- function( string )
{
    string <- paste0( string, collapse=" " )
    string <- gsub("___ ", "___", string, fixed=TRUE )
    form <- stats::as.formula(string)
    return(form)
}
