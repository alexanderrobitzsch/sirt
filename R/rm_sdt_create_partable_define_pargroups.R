## File Name: rm_sdt_create_partable_define_pargroups.R
## File Version: 0.09


rm_sdt_create_partable_define_pargroups <- function(partable, pg1, pg2)
{
    partable$pargroup <- 0
    # K <- max( partable$row, na.rm=TRUE)
    K <- max( partable$col, na.rm=TRUE)  # changed ARb 2019-07-21
    for (kk in 1:K){
        m1 <- max(partable$pargroup) + 1
        ind <- ( partable$type==pg1 ) & ( partable$col==kk)
        partable[ ind, "pargroup"] <- m1 * ( sum( partable[ind,"est"] ) > 0 )
    }
    m1 <- max(partable$pargroup) + 1
    ind <- ( partable$type==pg2 )
    partable[ ind, "pargroup"] <- m1 * ( sum( partable[ind,"est"] ) > 0 )
    return(partable)
}
