## File Name: xxirt_compute_prob_item_from_x.R
## File Version: 0.06

xxirt_compute_prob_item_from_x <- function(x, em_args, item_index=NULL)
{
    partable <- xxirt_partable_include_freeParameters( partable=em_args$partable,
                                        x=x[ em_args$parindex_items ] )
    probs_items <- xxirt_compute_itemprobs( item_list=em_args$item_list,
                                items=em_args$items, Theta=em_args$Theta,
                                ncat=em_args$ncat, partable=partable,
                                partable_index=em_args$partable_index,
                                item_index=item_index)
    return(probs_items)
}
