## File Name: MAR_normalized.R
## File Version: 0.02

MAR_normalized <- function(x)
{
    x <- stats::na.omit(x)
    MAR <- stats::median( abs( x - stats::median(x) ) )
    MAR <- MAR/0.6745
    return(MAR)
}
