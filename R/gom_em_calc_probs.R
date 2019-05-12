## File Name: gom_em_calc_probs.R
## File Version: 0.12



#--- gom calcprobs
gom_em_calc_probs <- function( lambda, theta.k, b=NULL, theta0.k=NULL )
{
    if ( ! is.null( b ) ){
        lambda <- stats::plogis( - b + matrix( theta0.k, length(b), length(theta0.k),
                    byrow=TRUE ) )
    }
    probs <- tcrossprod( lambda, theta.k)
    probsL <- array( 0, dim=c( nrow(lambda), 2, nrow(theta.k) ) )
    probsL[,2,] <- probs
    probsL[,1,] <- 1-probs
    res <- list("probs"=probs, "probsL"=probsL)
    return(res)
}

.gom.calcprobs <- gom_em_calc_probs
