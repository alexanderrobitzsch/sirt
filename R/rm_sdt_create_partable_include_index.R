## File Name: rm_sdt_create_partable_include_index.R
## File Version: 0.03


rm_sdt_create_partable_include_index <- function(partable)
{
    partable <- data.frame( index=1:nrow(partable), partable)
    return(partable)
}
