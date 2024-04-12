## File Name: xxirt_compute_itemprobs.R
## File Version: 0.215


# compute item probabilities
xxirt_compute_itemprobs <- function( item_list, items, Theta, ncat,
        partable, partable_index, item_index=NULL )
{
    TP <- nrow(Theta)
    maxK <- max(ncat)
    if ( is.null(item_index) ){
        I <- length(items)
        item_index <- 1L:I
    }
    I <- length(item_index)
    # compute item probabilities as a function of theta
    probs <- array( 0, dim=c(I,maxK,TP) )
    for (jj in 1L:I){
        ii <- item_index[jj]
        item_ii <- item_list[[ii]]
        par_ii <- partable[ partable_index[[ii]], 'value' ]
        ncat_ii <- ncat[ii]
        arg_ii <- list( par=par_ii, Theta=Theta, ncat=ncat_ii )
        probs_ii <- do.call( item_ii$P, arg_ii )
        probs[ jj, 1L:ncat_ii,] <- t(probs_ii)
    }
    return(probs)
}

