## File Name: linking_haberman_compute_median.R
## File Version: 0.02
## File Last Change: 2019-01-09

linking_haberman_compute_median <- function(x, w)
{
    res <- linking_haberman_remove_missings_vector(x=x,w=w)
    x <- res$x
    w <- res$w
    res <- stats::quantile(x=x, weights=w,    probs=.5, na.rm=TRUE, ties=TRUE)
    return(res)
}
