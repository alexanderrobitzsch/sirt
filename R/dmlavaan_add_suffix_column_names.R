## File Name: dmlavaan_add_suffix_column_names.R
## File Version: 0.02

dmlavaan_add_suffix_column_names <- function(x, suffix)
{
    colnames(x) <- paste0(colnames(x), suffix)
    return(x)
}
