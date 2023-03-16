## File Name: mgsem_add_list_entries.R
## File Version: 0.06


mgsem_add_list_entries <- function(list1, add_list, output_list, elements=NULL)
{
    if (is.null(elements)){
        elements <- names(add_list)
    }
    for (nn in elements){
        output_list[[nn]] <- list1[[nn]] + add_list[[nn]]
    }
    return(output_list)
}
