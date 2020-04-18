## File Name: invariance_alignment_process_parameters.R
## File Version: 0.03

invariance_alignment_process_parameters <- function(par.aligned, par)
{
    res <- sirt_matrix_names(x=par.aligned, extract_names=par)
    res[is.na(par)] <- NA
    return(res)
}
