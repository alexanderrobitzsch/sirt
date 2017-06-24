

##########################################################################
# prodterms function from mirt package
# This function is not exported and hence redefined in sirt
mirt_prodterms <- function (theta0, prodlist) 
{
    products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
    for (i in 1L:length(prodlist)) {
        tmp <- prodlist[[i]]
        for (j in 1L:length(tmp)){
			products[, i] <- products[, i] * theta0[, tmp[j]]
		}
    }
    ret <- cbind(theta0, products)
    return(ret)
}
