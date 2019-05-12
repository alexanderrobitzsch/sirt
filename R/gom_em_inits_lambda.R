## File Name: gom_em_inits_lambda.R
## File Version: 0.01

gom_em_inits_lambda <- function(I, K)
{
    lambda <- matrix( .75*seq( 1/(2*K), 1, 1/K), I, K, byrow=TRUE )
    return(lambda)
}
