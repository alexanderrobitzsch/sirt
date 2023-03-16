## File Name: xxirt_ic_compute_criteria.R
## File Version: 0.05
## File Last Change: 2022-01-26

xxirt_ic_compute_criteria <- function(ic, compute_np=TRUE)
{
    dev <- ic$deviance
    if (compute_np){
        ic$np <- ic$np.item + ic$np.Theta
    }
    ic$AIC <- dev + 2*ic$np
    log_n <- log(ic$n)
    ic$BIC <- dev + log_n * ic$np
    ic$CAIC <- dev + ( log_n + 1 )*ic$np
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )
    return(ic)
}
