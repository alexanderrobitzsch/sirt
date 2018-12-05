## File Name: tracemat.R
## File Version: 0.05

########################################
# trace of a matrix
tracemat <- function(A)
{
    res <- sum(diag(A))
    return(res)
}
########################################
