## File Name: gom_em_ic.R
## File Version: 0.06

gom_em_ic <- function(dev, dat2, I, K, TP, model, weights, lambda_partable=NULL)
{
    ic <- list( deviance=dev, n=sum(weights) )
    ic$np.item <- I*K
    if (!is.null(lambda_partable)){
        ic$np.item <- max(lambda_partable$par_index)
    }
    if (model=="GOMRasch"){
        ic$np.item <- I
    }
    # trait matrix
    ic$np.trait <- TP - 1
    if (model=="GOMnormal"){
        ic$np.trait <- 2*(K-1) + (K-1)*(K-2)/2
    }
    if (model=="GOMRasch"){
        ic$np.trait <- 3 + 1
    }
    ic$np <- ic$np.item + ic$np.trait
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
    # corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )
    #--- output
    return(ic)
}
