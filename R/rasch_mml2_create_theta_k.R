## File Name: rasch_mml2_create_theta_k.R
## File Version: 0.06
## File Last Change: 2021-05-02

rasch_mml2_create_theta_k <- function(Qmatrix, theta.k)
{
    if ( ( ! is.null( Qmatrix ) ) &  is.vector(theta.k) ){
        D <- ncol(Qmatrix)
        if ( D==2){ theta.k <- expand.grid( theta.k, theta.k ) }
        if ( D==3){ theta.k <- expand.grid( theta.k, theta.k, theta.k) }
        if ( D==4){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k) }
        if ( D==5){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k, theta.k) }
        if ( D==6){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k, theta.k,
                                    theta.k) }
        if ( D==7){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k, theta.k,
                                    theta.k, theta.k) }
        if ( D==8){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k, theta.k,
                                    theta.k, theta.k, theta.k) }
        if ( D==9){ theta.k <- expand.grid( theta.k, theta.k, theta.k, theta.k, theta.k,
                                    theta.k, theta.k, theta.k, theta.k) }
    }
    return(theta.k)
}
