## File Name: rm_sdt_create_partable_include_index.R
## File Version: 0.05
## File Last Change: 2018-12-30


rm_sdt_create_partable_include_index <- function(partable)
{
    partable <- data.frame( index=1:nrow(partable), partable)
    return(partable)
}
