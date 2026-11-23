## File Name: hh_linking_robust_weights_bisquare.R
## File Version: 0.01

hh_linking_robust_weights_bisquare <- function(x, k=1)
{
    res <- ifelse( abs(x) < k, (1-(x/k)^2)^2, 0 )
    return(res)
}
