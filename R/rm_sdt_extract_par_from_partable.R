## File Name: rm_sdt_extract_par_from_partable.R
## File Version: 0.01


rm_sdt_extract_par_from_partable <- function(partable)
{
    return( partable[ partable$est, "value" ] )
}
