## File Name: sirt_dmvnorm.R
## File Version: 0.113

sirt_dmvnorm <- function(x, mean=NULL, sigma=NULL, ... )
{
    TAM::require_namespace_msg("mvtnorm")
    if (is.matrix(x)){
        p <- ncol(x)
    } else {
        p <- 1
    }
    if (is.null(mean)){
        mean <- rep(0,p)
    }
    if (is.null(sigma)){
        sigma <- diag(p)
        if (p==1){
            sigma <- 1
        }
    }
    if (( p>1 ) | (is.matrix(x)) ){
        y <- mvtnorm::dmvnorm(x=x, mean=mean, sigma=sigma, ...)
    } else {
        y <- stats::dnorm(x=x, mean=mean, sd=sigma, ...)
    }
    return(y)
}
