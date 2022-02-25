## File Name: xxirt_ic.R
## File Version: 0.185


#-- information criteria xxirt
xxirt_ic <- function( dev, N, par_items, par_Theta, I, par_items_bounds )
{
    # Information criteria
    ic <- list( "deviance"=dev, "n"=N, "I"=I )
    ic$np.items <- sum(par_items_bounds$active)
    ic$np.Theta <- length(par_Theta)
    ic <- xxirt_ic_compute_criteria(ic=ic)
    return(ic)
}

