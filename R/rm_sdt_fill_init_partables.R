## File Name: rm_sdt_fill_init_partables.R
## File Version: 0.04
## File Last Change: 2018-12-30


rm_sdt_fill_init_partables <- function(partable_item, partable_rater,
    par_index, tau.item, a.item, c.rater, d.rater, type)
{
    partable_item <- rm_sdt_fill_init_partable( partable=partable_item,
                par_index=par_index, tau.item=tau.item, a.item=a.item, c.rater=c.rater,
                d.rater=d.rater, type='item' )
    partable_rater <- rm_sdt_fill_init_partable( partable=partable_rater,
                par_index=par_index, tau.item=tau.item, a.item=a.item, c.rater=c.rater,
                d.rater=d.rater, type='rater' )
    #--- output
    res <- list(partable_item=partable_item, partable_rater=partable_rater)
    return(res)
}
