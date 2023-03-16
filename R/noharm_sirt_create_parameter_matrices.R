## File Name: noharm_sirt_create_parameter_matrices.R
## File Version: 0.03
## File Last Change: 2019-01-05

noharm_sirt_create_parameter_matrices <- function(mat_label, parm_table, parm_index)
{
    nrow <- parm_index[[mat_label]]$nrow
    ncol <- parm_index[[mat_label]]$ncol
    mat <- matrix(0, nrow=nrow, ncol=ncol)
    x1 <- parm_table[ parm_index[[mat_label]]$row_parm_table, "est"]
    mat[ as.matrix(parm_index[[mat_label]]$entries) ] <- x1
    return(mat)
}
