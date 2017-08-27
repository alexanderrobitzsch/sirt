## File Name: diag2.R
## File Version: 0.02
## File Last Change: 2017-01-18 11:02:46

diag2 <- function( vec)
{
	if ( length(vec) > 1){
		res <- diag(vec)
	} else {
		res <- matrix(vec, nrow=1,ncol=1)
	}
	return(res)
}
