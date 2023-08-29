## File Name: dmlavaan_se_sandwich.R
## File Version: 0.115

dmlavaan_se_sandwich <- function(mod1, mod2, partable, label_parnames="parnames0",
        label_NPU="NPU", label_B="B", is_dmlavaan=TRUE)
{
    parnames1 <- mod1[[ label_parnames ]]
    NP1 <- mod1[[ label_NPU ]]
    parnames2 <- mod2[[ label_parnames ]]
    NP2 <- mod2[[ label_NPU ]]

    #-- compute standard error for difference
    scores <- cbind( mod1$scores, mod2$scores )
    N <- mod1$N
    A <- crossprod(scores)
    if (! is_dmlavaan){
        A <- A / N^4
    }
    B <- matrix(0, nrow=NP1+NP2, ncol=NP1+NP2)
    B[ 1:NP1, 1:NP1 ] <- mod1[[ label_B ]]
    B[ NP1 + 1:NP2, NP1 + 1:NP2 ] <- mod2[[ label_B ]]
    parnames0 <- c( paste0(parnames1, '_mod1'), paste0(parnames2, '_mod2') )
    res <- dmlavaan_sandwich_formula(A=A, B=B, parnames=parnames0)
    V <- res$V

    #-- now compute standard error for inference
    pars_sel <- which( ! is.na( partable$diff ) )

    pp <- pars_sel[1]
    for (pp in pars_sel){
        partable_pp <- partable[pp,]
        parname_pp <- partable_pp$parname
        ind_pp <- c( partable_pp$in_mod1, NP1 + partable_pp$in_mod2 )
        V_pp <- V[ ind_pp, ind_pp ]
        H_pp <- matrix( c(1,-1), nrow=1, ncol=2)
        DV_pp <- H_pp %*% V_pp %*% t(H_pp)
        partable$se_diff[pp] <- sqrt( DV_pp[1,1] )
    }
    partable$t <- partable$diff / partable$se_diff
    partable$p <- 2*stats::pnorm(-abs(partable$t))

    #--- output
    res <- list(partable=partable, V=V)
    return(res)
}
