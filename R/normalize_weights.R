## File Name: normalize_weights.R
## File Version: 0.01

normalize_weights <- function(x)
{
    N <- length(x)
    x <- N*x/sum(x)
    return(x)
}
