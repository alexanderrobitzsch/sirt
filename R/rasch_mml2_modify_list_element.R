## File Name: rasch_mml2_modify_list_element.R
## File Version: 0.01
## File Last Change: 2017-06-17 18:04:00

rasch_mml2_modify_list_element <- function( x, entry, value )
{
	x[[ entry ]] <- value
	return( x )
}
