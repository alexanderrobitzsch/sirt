
rasch_mml2_modify_list_element <- function( x, entry, value )
{
	x[[ entry ]] <- value
	return( x )
}