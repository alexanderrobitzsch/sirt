## File Name: rm_sdt_create_partable_include_fixed_item_category_parameters.R
## File Version: 0.11

rm_sdt_create_partable_include_fixed_item_category_parameters <- function(
    partable, tif, index_no_est_value )
{
    if ( ! is.null(tif) ){
        NT <- nrow(tif)
        for (tt in seq_len(NT) ){
            ind_tt <- which( ( partable$row==tif[tt,"item"] ) &
                                ( partable$col==tif[tt,"categ"] ) )
            partable[ ind_tt, "est"] <- FALSE
            partable[ ind_tt, "fixed"] <- TRUE
            partable[ ind_tt, "value"] <- tif[tt,"val"]
        }
        par_est <- unique( partable$parindex[ partable$est ] )
        partable$parindex <- match( partable$parindex, par_est )
        partable[ partable$fixed, "parindex" ] <- index_no_est_value
    }
    return(partable)
}
