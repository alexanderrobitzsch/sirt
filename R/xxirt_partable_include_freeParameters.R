## File Name: xxirt_partable_include_freeParameters.R
## File Version: 0.11
## File Last Change: 2022-01-26



xxirt_partable_include_freeParameters <- function( partable, x )
{
    vals <- x[ partable$parindex ]
    ind <- is.na(vals)
    vals[ ind ] <- partable$value[ ind ]
    partable$value <- vals
    return(partable)
}

