## File Name: linking_haberman_remove_missings_vector.R
## File Version: 0.01

linking_haberman_remove_missings_vector <- function(x,w)
{
    ind <- ! is.na(x)
    x <- x[ind]
    w <- w[ind]
    res <- list(x=x, w=w)
    return(res)
}
