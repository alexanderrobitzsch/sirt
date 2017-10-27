## File Name: linking_haberman_vcov_transformation.R
## File Version: 0.05

linking_haberman_vcov_transformation <- function( H1 , aj_vcov )
{
	aj_vcov <- H1 %*% aj_vcov %*% t(H1)
	aj_se <- c( sqrt( diag( aj_vcov ) ) )
	res <- list( vcov = aj_vcov , se = aj_se )
	return(res)
}
