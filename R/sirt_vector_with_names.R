## File Name: sirt_vector_with_names.R
## File Version: 0.06
## File Last Change: 2018-12-30

sirt_vector_with_names <- function(value, names)
{
    vec <- rep( value, length(names) )
    if ( ! is.numeric(names) ){
        names(vec) <- names
    }
    return(vec)
}
