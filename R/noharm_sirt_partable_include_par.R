## File Name: noharm_sirt_partable_include_par.R
## File Version: 0.02

noharm_sirt_partable_include_par <- function(par, parm_table)
{
    non_fixed <- attr(parm_table, 'non_fixed')
    include_index <- attr(parm_table, 'include_index')
    parm_table$est[ non_fixed ] <- par[ include_index ]
    return(parm_table)
}
