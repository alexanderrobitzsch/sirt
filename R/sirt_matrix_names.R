## File Name: sirt_matrix_names.R
## File Version: 0.07

sirt_matrix_names <- function(x, row_names=NULL, col_names=NULL, extract_names=NULL)
{
    if ( ! is.null(extract_names) ){
        row_names <- rownames(extract_names)
        col_names <- colnames(extract_names)
    }
    if ( ! is.null(row_names) ){
        rownames(x) <- row_names
    }
    if ( ! is.null(col_names) ){
        colnames(x) <- col_names
    }
    return(x)
}
