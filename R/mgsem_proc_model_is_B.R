## File Name: mgsem_proc_model_is_B.R
## File Version: 0.02


mgsem_proc_model_is_B <- function(model)
{
    H <- length(model)
    is_B <- 0
    for (hh in 1:H){
        is_B <- is_B + ( ! is.null( model[[hh]][["est"]][["B"]] ) )
    }
    is_B <- ( is_B > 0 )
    return(is_B)
}
