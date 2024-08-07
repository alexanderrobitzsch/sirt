## File Name: rm_facets_center_value_aggregate.R
## File Version: 0.071

rm_facets_center_value_aggregate <- function(x, index, value=0)
{
    y <- rep(NA, length(x) )
    index_unique <- unique(index)
    NI <- length(index_unique)
    for (ii in 1L:NI){
        index_ii <- index_unique[ii]
        ind_ii <- which( index==index_ii )
        y[ind_ii] <- rm_facets_center_value(x=x[ind_ii], value=value)
    }
    return(y)
}
