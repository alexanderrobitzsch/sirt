## File Name: meas_inv_compute_lavaan_parnames.R
## File Version: 0.01
## File Last Change: 2022-04-17


meas_inv_compute_lavaan_parnames <- function(object)
{
    pars <- paste0( object$lhs, object$op, object$rhs )
    return(pars)
}
