## File Name: bounds_parameters.R
## File Version: 0.08
## File Last Change: 2018-12-30


bounds_parameters <- function( pars, lower=NULL, upper=NULL)
{
    if ( ! is.null(lower)){
        pars <- ifelse( pars < lower, lower, pars )
    }
    if ( ! is.null(upper)){
        pars <- ifelse( pars > upper, upper, pars )
    }
    return(pars)
}
