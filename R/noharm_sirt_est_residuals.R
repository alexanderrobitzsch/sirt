## File Name: noharm_sirt_est_residuals.R
## File Version: 0.183


#**** estimate residuals
noharm_sirt_est_residuals <- function( Fval, Pval, Fpatt, Ppatt,
    I, D,  b0.jk, b1.jk, b2.jk, b3.jk, wgtm, pm, Psival, Psipatt, modesttype )
{
    # compute dj
    dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )

    # compute ej
    ej <- sqrt( 1 + dj^2 )
    ej.ek <- 1 / TAM::tam_outer( ej, ej )
    diag( ej.ek ) <- 0
    v0.jk <- b0.jk
    v1.jk <- b1.jk * ej.ek
    v2.jk <- b2.jk * ej.ek^2
    v3.jk <- b3.jk * ej.ek^3
    # compute gamma.jk=f_j' P f_k
    gamma.jk <- Fval %*% Pval %*% t(Fval) + Psival
    # gamma.jk <- gamma.jk * ej.ek
    # compute p_d ' f_k
    pd.fk <- Fval %*% Pval
    Fval_old <- Fval
    Pval_old <- Pval
    residuals <- ( wgtm * ( pm - v0.jk - v1.jk*gamma.jk -
                                    v2.jk*gamma.jk^2 - v3.jk*gamma.jk^3 ) )

    #* fit statistics
    rmsr <- sqrt( sum( residuals^2 * wgtm ) / sum(wgtm) )
    RMW <- residuals * wgtm
    PMW <- pm * wgtm
    tanaka <- 1 - sum( diag( RMW %*% RMW ) ) /  sum(diag( PMW %*% PMW ))
    if (modesttype==2){
        tanaka <- rmsr <- NA
    }

    #--- output
    res <- list(residuals=residuals, rmsr=rmsr, tanaka=tanaka)
    return(res)
}
