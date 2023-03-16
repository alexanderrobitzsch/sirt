## File Name: weighted_rowMeans.R
## File Version: 0.09
## File Last Change: 2018-12-30

weighted_rowMeans <- function( mat, wgt=NULL)
{
    wgt <- weighted_stats_extend_wgt( wgt=wgt, mat=mat )
    mat1 <- rowSums( mat * wgt, na.rm=TRUE)
    mat2 <- rowSums( wgt, na.rm=TRUE)
    mat1 <- mat1 / mat2
    return(mat1)
}
