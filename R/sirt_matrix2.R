
sirt_matrix2 <- function(x , nrow)
{
	matr <- matrix( x , nrow=nrow, ncol=length(x) , byrow=TRUE )
	return(matr)
}