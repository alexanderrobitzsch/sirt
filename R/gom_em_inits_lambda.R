## File Name: gom_em_inits_lambda.R
## File Version: 0.02

gom_em_inits_lambda <- function(I, K, lambda.inits=NULL)
{
    if (is.null(lambda.inits)){
        lambda <- matrix( .75*seq( 1/(2*K), 1, 1/K), I, K, byrow=TRUE )
    } else {
        lambda <- lambda.inits
    }
    return(lambda)
}
