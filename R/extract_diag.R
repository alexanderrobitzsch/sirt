## File Name: extract_diag.R
## File Version: 0.02

extract_diag <- function(X, indices=NULL)
{
    if (is.null(indices)){
        indices <- 1L:nrow(X)
    }
    d <- diag( X[indices,indices] )
    return(d)
}
