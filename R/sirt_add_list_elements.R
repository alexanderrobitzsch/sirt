## File Name: sirt_add_list_elements.R
## File Version: 0.04

sirt_add_list_elements <- function(res, res2)
{
    NR <- length(res2)
    for (rr in 1L:NR){
        res[[ names(res2)[[rr]] ]] <- res2[[rr]]
    }
    return(res)
}
