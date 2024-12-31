## File Name: sqrt_diag.R
## File Version: 0.03

sqrt_diag <- function(x, names=NULL)
{
    v1 <- diag(x)
    v1 <- v1*(v1>0)
    res <- sqrt(v1)
    names(res) <- colnames(x)
    if (!is.null(names)){
        names(res) <- names
    }
    return(res)
}
