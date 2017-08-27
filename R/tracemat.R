## File Name: tracemat.R
## File Version: 0.03
## File Last Change: 2017-01-18 11:02:55

########################################
# trace of a matrix
tracemat <- function(A)
{ 
    res <- sum( diag(A) ) 
	return(res)
} 
########################################
