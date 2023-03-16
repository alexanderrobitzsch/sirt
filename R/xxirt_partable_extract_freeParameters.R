## File Name: xxirt_partable_extract_freeParameters.R
## File Version: 0.11
## File Last Change: 2022-01-26


xxirt_partable_extract_freeParameters <- function( partable )
{
    partable <- partable[ partable$est, ]
    partable <- partable[ order(partable$parindex), ]
    partable <- partable[ ! duplicated( partable$parindex), ]
    x <- partable$value
    names(x) <- partable$parlabel
    return(x)
}

