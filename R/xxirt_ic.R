## File Name: xxirt_ic.R
## File Version: 0.195


#-- information criteria xxirt
xxirt_ic <- function( dev, N, par_items, par_Theta, I, par_items_bounds, np_item=NULL,
            np_theta=NULL, estimator="ML")
{
    # Information criteria
    ic <- list( deviance=dev, n=N, I=I )
    ic$np.items <- sum(par_items_bounds$active)
    if ( ! is.null(np_item) ){
        ic$np.items <- np_item
    }
    NPT <- length(par_Theta)
    if (!is.null(np_theta)){
        NPT <- np_theta
    }
    ic$np.Theta <- NPT
    ic <- xxirt_ic_compute_criteria(ic=ic, estimator=estimator)
    return(ic)
}

