## File Name: sirt_matrix2.R
## File Version: 0.01
## File Last Change: 2017-06-17 17:00:03

sirt_matrix2 <- function(x , nrow)
{
	matr <- matrix( x , nrow=nrow, ncol=length(x) , byrow=TRUE )
	return(matr)
}
