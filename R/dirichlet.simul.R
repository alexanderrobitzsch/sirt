## File Name: dirichlet.simul.R
## File Version: 0.131


#-- simulate from a Dirichlet distribution
dirichlet.simul <- function( alpha )
{
    N <- nrow(alpha)
    K <- ncol(alpha)
    ygamma <- 0*alpha
    for (ii in 1L:K){
        ygamma[,ii] <- stats::rgamma( n=N, shape=alpha[,ii] )
    }
    x <- ygamma / rowSums(ygamma)
    return(x)
}
