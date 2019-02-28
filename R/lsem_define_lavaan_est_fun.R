## File Name: lsem_define_lavaan_est_fun.R
## File Version: 0.03

lsem_define_lavaan_est_fun <- function(lavaan_fct)
{
    if (lavaan_fct=="sem"){
        lavaan_est_fun <- lavaan::sem
    }
    if (lavaan_fct=="lavaan"){
        lavaan_est_fun <- lavaan::lavaan
    }
    if (lavaan_fct=="cfa"){
        lavaan_est_fun <- lavaan::cfa
    }
    return(lavaan_est_fun)
}
