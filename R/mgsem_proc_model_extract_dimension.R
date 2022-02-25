## File Name: mgsem_proc_model_extract_dimension.R
## File Version: 0.09


mgsem_proc_model_extract_dimension <- function(model, entry="est", type, nrow=TRUE)
{
    G <- length(model)
    gg <- 1
    while (gg<=G){
        mat <- model[[gg]][[entry]][[type]]
        if (!is.null(mat)){
            val <- ifelse(nrow, nrow(mat), ncol(mat) )
            gg <- G+1
        }
        gg <- gg + 1
    }
    return(val)
}

