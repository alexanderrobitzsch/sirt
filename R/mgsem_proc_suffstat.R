## File Name: mgsem_proc_suffstat.R
## File Version: 0.097


mgsem_proc_suffstat <- function(suffstat)
{
    G <- length(suffstat)
    N_group <- rep(0,G)
    for (gg in 1:G){
        suffstat_gg <- suffstat[[gg]]
        N_group[gg] <- suffstat_gg$N
        # rename entries
        ns_gg <- names(suffstat[[gg]])
        if( ! ( "M" %in% ns_gg ) ){
            ind <- min( which( unlist( lapply(suffstat_gg, is.vector) ) ))
            ns_gg[ind] <- "M"
        }
        if( ! ( "S" %in% ns_gg ) ){
            ind <- min( which( unlist( lapply(suffstat_gg, is.matrix) ) ))
            ns_gg[ind] <- "S"
        }
        names(suffstat_gg) <- ns_gg
        I <- length(suffstat_gg$M)

        #*** weights for moment-based estimation
        if ( ! ( "weights_M" %in% ns_gg ) ){
            suffstat_gg$weights_M <- 1+0*suffstat_gg$M
        }
        if ( ! ( "weights_S" %in% ns_gg ) ){
            suffstat_gg$weights_S <- 1+0*suffstat_gg$S
        }
        suffstat_gg$vech_weights_S <- mgsem_vech(x=suffstat_gg$weights_S)
        suffstat[[gg]] <- suffstat_gg
    }  # end gg
    N <- sum(N_group)

    #--- output
    res <- list(G=G, N=N, I=I, N_group=N_group, suffstat=suffstat)
    return(res)
}
