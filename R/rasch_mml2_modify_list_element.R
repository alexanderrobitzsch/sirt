## File Name: rasch_mml2_modify_list_element.R
## File Version: 0.03
## File Last Change: 2018-12-30

rasch_mml2_modify_list_element <- function( x, entry, value )
{
    x[[ entry ]] <- value
    return( x )
}
