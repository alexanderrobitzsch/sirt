## File Name: rm_sdt_extract_par_from_partable.R
## File Version: 0.03
## File Last Change: 2018-12-30


rm_sdt_extract_par_from_partable <- function(partable)
{
    return( partable[ partable$est, "value" ] )
}
