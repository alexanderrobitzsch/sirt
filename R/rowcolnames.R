## File Name: rowcolnames.R
## File Version: 0.02

rowcolnames <- function(x, names)
{
    rownames(x) <- colnames(x) <- names
    return(x)
}
