## File Name: rm_sdt_create_parm_index_rater.R
## File Version: 0.04


rm_sdt_create_parm_index_rater <- function( est.rater, ND, item.index,
    rater.index)
{
    g1 <- NULL
    if ( est.rater=="a"){ g1 <- seq_len(ND) }
    if ( est.rater=="r"){ g1 <- rater.index }
    if ( est.rater=="i"){ g1 <- item.index }
    if ( est.rater=="e"){ g1 <- rep(1,ND) }
    if ( est.rater=="n"){ g1 <- rep(-999,ND) }
    return(g1)
}
