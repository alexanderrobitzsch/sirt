## File Name: sirt_pem_extract_parameters.R
## File Version: 0.01
## File Last Change: 2017-10-03 11:28:14

sirt_pem_extract_parameters <- function( parm, parmgroup, pem_parameter_index )
{
	info <- pem_parameter_index[[ parmgroup ]]
	x_dim <- info$dim
	x <- parm[ info$index ]
	x <- sirt_pem_adjust_dimension(x=x, x_dim=x_dim )
	return(x)
}
