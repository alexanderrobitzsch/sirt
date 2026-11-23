## File Name: hh_linking_apply_robust_weights.R
## File Version: 0.087

hh_linking_apply_robust_weights <- function(residuals, rob_weights="none", V=NULL,
        eps=0.001, c=NULL)
{
    k <- NULL
    robw <- rep(1, length(residuals))
    if (rob_weights %in% c('bisquare1','bisquare2')){
        k <- hh_linking_robust_weights_bisquare_tuning_constant(x=residuals,
                    rob_weights=rob_weights, V=V, c=c)
        c <- attr(k, 'c')
        robw <- hh_linking_robust_weights_bisquare(x=residuals, k=k)
    }
    if (substring(rob_weights,1,1)=='L'){
        p <- as.numeric(paste(substring(rob_weights,2)))
        robw <- hh_linking_robust_weights_Lp(x=residuals, p=p, eps=eps)
    }
    robw <- normalize_weights(x=robw)
    attr(robw,'k') <- k

    #- output
    return(robw)
}
