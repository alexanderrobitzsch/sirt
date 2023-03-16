## File Name: dmlavaan_add_suffix_column_names.R
## File Version: 0.02
## File Last Change: 2023-03-10

dmlavaan_add_suffix_column_names <- function(x, suffix)
{
    colnames(x) <- paste0(colnames(x), suffix)
    return(x)
}
