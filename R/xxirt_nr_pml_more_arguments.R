## File Name: xxirt_nr_pml_more_arguments.R
## File Version: 0.129

xxirt_nr_pml_more_arguments <- function(em_args, grad_fun)
{
    parindex <- em_args$partable$parindex
    parindex1 <- setdiff(parindex,0)
    if( length(parindex1)>length(unique(parindex1))){
        grad_fun <- NULL
        no_itempar_constraints <- FALSE
    } else {
        no_itempar_constraints <- TRUE
    }
    parindex_items <- em_args$parindex_items
    parindex_Theta <- em_args$parindex_Theta
    em_args$pml_args$in_par_Theta <- sapply(1L:em_args$NP, FUN=function(pp){
                                                pp %in% parindex_Theta },
                                                simplify=TRUE )
    #- index handling
    I <- em_args$I
    NPI <- em_args$NPI
    NP <- em_args$NP
    NI2 <- em_args$pml_args$NI2

    index_freq2 <- index_freq1 <- as.list(1L:NP)
    W2_long <- em_args$pml_args$W2_long

    for (pp in 1L:NP){
        if (pp<=NPI){
            item_pp <- em_args$item_index[[pp]]
            index_freq1[[pp]] <- item_pp - 1
            v1 <- sort( union( which( W2_long[,1] %in% item_pp ),
                                which( W2_long[,2] %in% item_pp ) ) )
            index_freq2[[pp]] <- v1 - 1
        } else {
            index_freq1[[pp]] <- (1L:I) - 1
            index_freq2[[pp]] <- (1L:NI2) - 1
        }
    }
    em_args$pml_args$index_freq1 <- index_freq1
    em_args$pml_args$index_freq2 <- index_freq2
    em_args$pml_args$no_itempar_constraints <- no_itempar_constraints

    #-- output
    res <- list(em_args=em_args, grad_fun=grad_fun)
    return(res)
}
