## File Name: rm.facets_alg.R
## File Version: 4.28

#######################################################
# parameters expanded dataset
rm.facets.itempar.expanded <- function( b.item, b.rater, Qmatrix, tau.item,
        VV, K, I, TP, a.item, a.rater, item.index, rater.index,
        theta.k, RR )
{
    b <- tau.item[ item.index, ]
    b0 <- ( matrix( b.rater, nrow=RR, ncol=K) )[ rater.index, ] *     Qmatrix[ item.index,]
    b <- b + b0
    # a parameter
    a <- a.item[ item.index ] * a.rater[ rater.index ]
    res <- list("a"=a, "b"=b )
     return(res)
    }
#########################################################



########################################################
# calculation of probabilities in the partial credit model
.rm.pcm.calcprobs <- function( a, b, Qmatrix, theta.k, I, K, TP )
{
    probs <- array( 0, dim=c(I,K+1,TP) )   # categories 0, ..., K
    for (kk in 1:K){
        l0 <- matrix( - b[,kk], nrow=I,ncol=TP)
        l0 <- l0 + outer( a * Qmatrix[, kk], theta.k )
        probs[,kk+1,] <- l0
    }
    probs <- exp( probs )
    probs1 <- probs[,1,]
    for (kk in 2:(K+1)){
        probs1 <- probs1 + probs[,kk,]
    }
    for (kk in 1:(K+1)){
        probs[,kk,] <- probs[,kk,] / probs1
    }
   return(probs)
}

######################
sumtau <- function(tau.item)
{
    K <- ncol(tau.item)
    matr <- tau.item
    for (kk in 2:K){
        matr[,kk] <- rowSums( tau.item[,1:kk] )
    }
    return(matr)
}

#############################################################################
# calculation of probabilities in the facet model
.rm.facets.calcprobs <- function( b.item, b.rater, Qmatrix, tau.item,
        VV, K, I, TP, a.item, a.rater, item.index, rater.index,
        theta.k, RR )
{
    b <- tau.item[ item.index, ]
    b0 <- ( matrix( b.rater, nrow=RR, ncol=K) )[ rater.index, ] * Qmatrix[ item.index,]
    b <- b + b0
    # a parameter
    a <- a.item[ item.index ] * a.rater[ rater.index ]
    res <- .rm.pcm.calcprobs( a, b, Qmatrix=Qmatrix[item.index,], theta.k, I, K, TP )
    return(res)
}


