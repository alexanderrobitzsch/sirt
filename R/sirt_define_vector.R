## File Name: sirt_define_vector.R
## File Version: 0.01
## File Last Change: 2019-04-26

sirt_define_vector <- function( value, names)
{
    NN <- length(names)
    res <- rep(value, NN)
    names(res) <- names
    return(res)
}
