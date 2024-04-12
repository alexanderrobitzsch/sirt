## File Name: sirt_matrix_lower_to_upper.R
## File Version: 0.03

sirt_matrix_lower_to_upper <- function(x)
{
    N <- ncol(x)
    for (ii in 1L:N){
        for (jj in ii:N){
            x[ii,jj] <- x[jj,ii]
        }
    }
    return(x)
}
