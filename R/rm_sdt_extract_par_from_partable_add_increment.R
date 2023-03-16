## File Name: rm_sdt_extract_par_from_partable_add_increment.R
## File Version: 0.07
## File Last Change: 2018-12-30


rm_sdt_extract_par_from_partable_add_increment <- function(partable,
        pargroup, increment )
{
    partable <- partable[ partable$est, ]
    res <- partable[ partable$est, "value" ]
    res <- res + increment * ( partable$pargroup==pargroup )
    return(res)
}
