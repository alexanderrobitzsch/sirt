## File Name: sirt_vector_with_names.R
## File Version: 0.02

sirt_vector_with_names <- function(value, names)
{
    vec <- rep( value, length(names) )
    names(vec) <- names
    return(vec)
}
