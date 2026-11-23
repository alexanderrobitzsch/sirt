## File Name: hh_linking_robust_weights_bisquare_tuning_constant.R
## File Version: 0.087

hh_linking_robust_weights_bisquare_tuning_constant <- function(x, rob_weights, V=NULL,
        c=NULL)
{
    if (rob_weights=='bisquare1'){
        if (is.null(c)){
            c <- 1.96
        }
        k <- c*sqrt(diag(V))
    }
    if (rob_weights=='bisquare2'){
        if (is.null(c)){
            c <- 3
        }
        k <- c*stats::mad(x) # 4.685 is also a frequent choice
    }
    attr(k, 'c') <- c
    return(k)
}
