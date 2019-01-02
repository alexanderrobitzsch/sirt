## File Name: weighted_rowSums.R
## File Version: 0.08


weighted_rowSums <- function( mat, wgt=NULL)
{
    wgt <- weighted_stats_extend_wgt( wgt=wgt, mat=mat )
    mat1 <- rowSums( mat * wgt, na.rm=TRUE)
    return(mat1)
}
