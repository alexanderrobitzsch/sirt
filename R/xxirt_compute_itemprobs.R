## File Name: xxirt_compute_itemprobs.R
## File Version: 0.226


# compute item probabilities
xxirt_compute_itemprobs <- function( item_list, items, Theta, ncat,
        partable, partable_index, item_index=NULL )
{
    person_covariates <- attr(partable, 'person_covariates')
    TP <- nrow(Theta)
    maxK <- max(ncat)
    if ( is.null(item_index) ){
        I <- length(items)
        item_index <- 1L:I
    }
    I <- length(item_index)
    # compute item probabilities as a function of theta
    if (person_covariates){
        N <- nrow( (item_list[[1]])$X )
        dim_probs <- c(I,maxK,TP,N)
    } else {
        dim_probs <- c(I,maxK,TP)
    }
    probs <- array( 0, dim=dim_probs )
    for (jj in 1L:I){
        ii <- item_index[jj]
        item_ii <- item_list[[ii]]
        par_ii <- partable[ partable_index[[ii]], 'value' ]
        ncat_ii <- ncat[ii]
        arg_ii <- list( par=par_ii, Theta=Theta, ncat=ncat_ii )
        if (person_covariates){
            arg_ii$X <- item_ii$X
        }
        probs_ii <- do.call( what=item_ii$P, args=arg_ii )
        if (!person_covariates){
            probs[ jj, 1L:ncat_ii,] <- t(probs_ii)
        } else {
            probs[jj,,,] <- aperm(probs_ii, perm=c(2,1,3) )
        }
    }  # end jj
    return(probs)
}

