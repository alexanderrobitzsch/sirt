## File Name: rm_facets_itempar_expanded.R
## File Version: 0.06



#######################################################
# parameters expanded dataset
rm_facets_itempar_expanded <- function( b.item, b.rater, Qmatrix, tau.item,
        VV, K, I, TP, a.item, a.rater, item.index, rater.index,
        theta.k, RR )
{
    b <- tau.item[ item.index, ]
    b0 <- ( matrix( b.rater, nrow=RR, ncol=K) )[ rater.index, ] * Qmatrix[ item.index,]
    b <- b + b0
    # a parameter
    a <- a.item[ item.index ] * a.rater[ rater.index ]
    res <- list(a=a, b=b )
    return(res)
}
#########################################################
