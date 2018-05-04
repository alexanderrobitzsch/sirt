## File Name: sirt_pem_extract_dimension.R
## File Version: 0.01

sirt_pem_extract_dimension <- function(x)
{
    x_dim <- NULL
    if ( is.vector(x) ){
        x_dim <- length(x)
    }
    if ( is.matrix(x) | is.array(x) | is.data.frame(x) ){
        x_dim <- dim(x)
    }
    return(x_dim)
}
