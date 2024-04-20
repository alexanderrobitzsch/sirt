## File Name: noharm_sirt_partable_extract_par.R
## File Version: 0.03

noharm_sirt_partable_extract_par <- function(parm_table, col="est")
{
    extract_index <- attr(parm_table, 'extract_index')
    x <- parm_table[ extract_index, col]
    return(x)
}
