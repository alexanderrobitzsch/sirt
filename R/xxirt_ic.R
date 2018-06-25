## File Name: xxirt_ic.R
## File Version: 0.14


###################################################################
xxirt_ic <- function( dev, N, par_items, par_Theta, I,
       par_items_bounds )
{
    # Information criteria
    ic <- list( "deviance"=dev, "n"=N, "I"=I )
    # ic$np.item <- length(par_items)
    ic$np.items <- sum(par_items_bounds$active)
    ic$np.Theta <- length(par_Theta)
    ic$np <- ic$np.item + ic$np.Theta
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    log_n <- log(ic$n)
    ic$BIC <- dev + log_n * ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( log_n + 1 )*ic$np
    # corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )
    return(ic)
}
###################################################################
