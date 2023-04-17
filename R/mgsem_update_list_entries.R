## File Name: mgsem_update_list_entries.R
## File Version: 0.041


mgsem_update_list_entries <- function(add_list, output_list, elements=NULL)
{
    if (is.null(elements)){
        elements <- names(add_list)
    }
    for (nn in elements){
        output_list[[nn]] <- add_list[[nn]]
    }
    return(output_list)
}
