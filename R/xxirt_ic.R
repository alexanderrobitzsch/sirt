## File Name: xxirt_ic.R
## File Version: 0.191


#-- information criteria xxirt
xxirt_ic <- function( dev, N, par_items, par_Theta, I, par_items_bounds, np_item=NULL )
{
    # Information criteria
    ic <- list( deviance=dev, n=N, I=I )
    ic$np.items <- sum(par_items_bounds$active)
    if ( ! is.null(np_item) ){
        ic$np.items <- np_item
    }
    ic$np.Theta <- length(par_Theta)
    ic <- xxirt_ic_compute_criteria(ic=ic)
    return(ic)
}

