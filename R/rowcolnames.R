## File Name: rowcolnames.R
## File Version: 0.02
## File Last Change: 2023-03-11

rowcolnames <- function(x, names)
{
    rownames(x) <- colnames(x) <- names
    return(x)
}
