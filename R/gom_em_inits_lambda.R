## File Name: gom_em_inits_lambda.R
## File Version: 0.05

gom_em_inits_lambda <- function(I, K, lambda.inits=NULL, lambda_partable=NULL)
{
    if (is.null(lambda.inits)){
        lambda <- matrix( seq( 1/(2*K), 1, 1/K), I, K, byrow=TRUE )
    } else {
        lambda <- lambda.inits
    }
    a1 <- stats::aggregate( as.vector(lambda), list(lambda_partable$par_index), mean )
    lambda <- matrix( a1[lambda_partable$par_index,2], nrow=I, ncol=K)
    return(lambda)
}
