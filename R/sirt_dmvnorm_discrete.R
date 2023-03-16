## File Name: sirt_dmvnorm_discrete.R
## File Version: 0.07

sirt_dmvnorm_discrete <- function(x, mean=NULL, sigma=NULL, as_matrix=FALSE,
    eps=0, ...)
{
    y <- sirt_dmvnorm(x=x, mean=mean, sigma=sigma, ... )
    if (as_matrix){
        y <- matrix(y, nrow=length(y), ncol=1)
    }
    y <- y + eps
    y <- y / sum(y)
    return(y)
}
