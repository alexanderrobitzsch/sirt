## File Name: sirt_define_vector.R
## File Version: 0.01

sirt_define_vector <- function( value, names)
{
    NN <- length(names)
    res <- rep(value, NN)
    names(res) <- names
    return(res)
}
