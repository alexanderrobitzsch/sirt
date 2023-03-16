## File Name: tracemat.R
## File Version: 0.07
## File Last Change: 2018-12-30

########################################
# trace of a matrix
tracemat <- function(A)
{
    res <- sum(diag(A))
    return(res)
}
########################################
