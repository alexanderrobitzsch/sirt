## File Name: invariance_alignment_proc_labels.R
## File Version: 0.06

invariance_alignment_proc_labels <- function(x)
{
    G <- nrow(x)
    I <- ncol(x)
    if (is.null(colnames(x))){
        colnames(x) <- paste0('I', 1L:I)
    }
    if (is.null(rownames(x))){
        rownames(x) <- paste0('G', 1L:G)
    }
    return(x)
}
