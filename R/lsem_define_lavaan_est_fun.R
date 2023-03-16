## File Name: lsem_define_lavaan_est_fun.R
## File Version: 0.101

lsem_define_lavaan_est_fun <- function(lavaan_fct)
{
    TAM::require_namespace_msg('lavaan')
    if (lavaan_fct=='sem'){
        lavaan_est_fun <- lavaan::sem
    }
    if (lavaan_fct=='lavaan'){
        lavaan_est_fun <- lavaan::lavaan
    }
    if (lavaan_fct=='cfa'){
        lavaan_est_fun <- lavaan::cfa
    }
    if (lavaan_fct=='growth'){
        lavaan_est_fun <- lavaan::growth
    }
    return(lavaan_est_fun)
}
