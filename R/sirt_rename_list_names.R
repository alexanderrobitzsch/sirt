## File Name: sirt_rename_list_names.R
## File Version: 0.03


sirt_rename_list_names <- function(x, old, new)
{
    ind1 <- which( names(x)==old )
    if (length(ind1)>0){
        names(x)[ind1] <- new
    }
    return(x)
}
