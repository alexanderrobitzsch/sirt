## File Name: noharm_sirt_compute_final_constants.R
## File Version: 0.05
## File Last Change: 2018-12-30

noharm_sirt_compute_final_constants <- function(Fval, Pval, betaj,
    modesttype)
{
    #---- calculate final constants
    if (modesttype==2){
        dj2 <- diag( Fval %*% Pval %*% t(Fval) )
        Fval <- Fval / sqrt( 1 - dj2 )
    }
    # recalculation of f0 coefficient (final constant)
    dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )
    ej <- sqrt( 1 + dj^2 )
    f0 <- - betaj * ej
    # uniquenesses
    uqn <- 1 - ( dj^2 / ( 1 + dj^2 ) )
    # standardized loadings
    loadingsF <- Fval / ej

    #--- output
    res <- list(f0=f0, uqn=uqn, loadingsF=loadingsF)
    return(res)
}
