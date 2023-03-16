## File Name: sirt_add_names.R
## File Version: 0.02

sirt_add_names <- function(x, names)
{
    if (! is.null(names)){
        if (is.vector(x)){
            names(x) <- names
        }
        if (is.matrix(x)){
            rownames(x) <- names
            colnames(x) <- names
        }
    }
    return(x)
}
