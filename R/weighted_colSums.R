## File Name: weighted_colSums.R
## File Version: 0.04


weighted_colSums <- function( mat , wgt=NULL)
{
    wgt <- weighted_stats_extend_wgt( wgt=wgt , mat=mat )
    mat1 <- colSums( mat * wgt , na.rm=TRUE) 
    return(mat1)
}
