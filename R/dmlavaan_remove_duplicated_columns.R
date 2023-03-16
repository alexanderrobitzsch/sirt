## File Name: dmlavaan_remove_duplicated_columns.R
## File Version: 0.03

dmlavaan_remove_duplicated_columns <- function(x)
{
    x <- x[, ! duplicated(colnames(x)) ]
    return(x)
}
