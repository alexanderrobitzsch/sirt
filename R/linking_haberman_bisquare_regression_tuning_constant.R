## File Name: linking_haberman_bisquare_regression_tuning_constant.R
## File Version: 0.01

linking_haberman_bisquare_regression_tuning_constant <- function(x)
{
    x <- stats::na.omit(x)
    MAR <- stats::median( abs( x - stats::median(x) ) )
    k <- 4.685*MAR/0.6745
    return(k)
}
