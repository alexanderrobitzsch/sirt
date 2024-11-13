## File Name: linking_haberman_lq_pw_le.R
## File Version: 0.132

linking_haberman_lq_pw_le <- function(des, res_optim, vcov_list=NULL, symm_hess=FALSE)
{
    requireNamespace('MASS')
    ind_studies <- des$ind_studies
    ind_items <- des$ind_items
    design <- des$design
    I <- des$I
    G <- des$G
    par_delta <- c( exp(res_optim$slopes$coefficients),
                        res_optim$intercepts$coefficients )
    parnames <- c( paste0('sig',2L:G), paste0('mu',2L:G) )
    names(par_delta) <- parnames
    par_gamma <- rep(NA, 2*I*G)
    for (gg in 1L:G){
        v1 <- paste0( rep( paste0('I', 1L:I), each=2), '_', rep(c('a','b'), I), '_G', gg )
        names(par_gamma)[ 2*I*(gg-1)+seq(1,2*I) ] <- v1
    }
    NIP <- nrow(itempars)
    for (hh in 1L:NIP){
        ind <- 2*I*(ind_studies[hh]-1)+ 2*(ind_items[hh]-1)+1
        par_gamma[ ind ] <- itempars[hh,3]
        par_gamma[ ind+1 ] <- itempars[hh,4]
    }

    #** arrange variance matrix of item parameters
    Vgamma <- linking_haberman_lq_pw_le_arrange_Vgamma( vcov_list=vcov_list,
                    par_gamma=par_gamma, I=I, G=G, ind_items=ind_items,
                    ind_studies=ind_studies )

    #*** compute derivatives
    args1 <- list(par_delta=par_delta, par_gamma=par_gamma, des=des)

    # H_delta
    grad_delta <- do.call(what=linking_haberman_lq_pw_le_grad, args=args1)
    # H_delta_delta
    hess_delta_delta <- do.call(what=linking_haberman_lq_pw_le_hess_delta, args=args1)
    if (symm_hess){
        hess_delta_delta <- ( hess_delta_delta + t(hess_delta_delta) ) / 2
    }
    HDD <- MASS::ginv(X=hess_delta_delta)

    # H_delta_gamma
    hess_delta_gamma <- do.call(what=linking_haberman_lq_pw_le_hess_gamma, args=args1)

    #-- compute standard errors
    pn2 <- list(parnames, parnames)
    A <- HDD %*% hess_delta_gamma
    V_SE <- A %*% Vgamma %*% t(A)
    dimnames(V_SE) <- pn2
    SE <- sqrt_diag_positive(x=V_SE)

    #-- compute linking errors
    B <- crossprod( grad_delta )
    V_LE <- I/(I-1) * HDD %*% B %*% t(HDD)
    dimnames(V_LE) <- pn2
    LE <- sqrt_diag_positive(x=V_LE)

    #-- compute bias-corrected linking errors
    V_LEbc <- 0*V_LE
    for (ii in 1L:I){
        ind_ii <- rep(2*I*((1L:G)-1),each=2) + rep( 2*(ii-1)+1L:2, G)
        V_gammai <- Vgamma[ind_ii, ind_ii]
        H_delta_gammai <- hess_delta_gamma[,ind_ii]
        V_LEbc <- V_LEbc + H_delta_gammai %*% V_gammai %*% t(H_delta_gammai)
    }
    V_LEbc <- I/(I-1) * HDD %*% V_LEbc %*% t(HDD)
    dimnames(V_LEbc) <- pn2
    V_LEbc <- V_LE - V_LEbc
    LEbc <- sqrt_diag_positive(x=V_LEbc)

    #-- compute total errors
    V_TE <- V_SE + V_LE
    dimnames(V_TE) <- pn2
    TE <- sqrt_diag_positive(x=V_TE)
    TEbc <- sqrt( SE^2 + LEbc^2 )

    #-- output
    res <- list(V_SE=V_SE, SE=SE, V_LE=V_LE, LE=LE, V_TE=V_TE, TE=TE,
                    V_LEbc=V_LEbc, LEbc=LEbc, TEbc=TEbc,
                    grad_delta=grad_delta, hess_delta_delta=hess_delta_delta,
                    hess_delta_gamma=hess_delta_gamma, I=I, G=G)
    return(res)
}
