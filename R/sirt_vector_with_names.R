## File Name: sirt_vector_with_names.R
## File Version: 0.01
## File Last Change: 2017-09-20 10:18:52

sirt_vector_with_names <- function(value, names)
{
	vec <- rep( value , length(names) )
	names(vec) <- names
	return(vec)
}
