## File Name: rm_sdt_fill_par_to_partable.R
## File Version: 0.18
## File Last Change: 2018-12-30


rm_sdt_fill_par_to_partable <- function( par_index, partable, parm0, type)
{

    not_fixed <- ! partable$fixed
    partable_parindex <- partable$parindex[ not_fixed ]
    partable[ not_fixed, "value"] <- parm0[ partable_parindex ]
    #--- fill objects
    parm_list <- NULL
    if (type=="item"){
        pars <- c("tau.item", "a.item")
    } else {
        pars <- c("c.rater", "d.rater")
    }
    for (pp in pars){
        index_pp <- par_index[[pp]]
        n_pp <- index_pp[ length(index_pp) ]
        n_row <- partable[n_pp, "row"]
        n_col <- partable[n_pp, "col"]
        if ( ! is.na(n_col) ){
            x <- matrix( partable[ index_pp, "value"], nrow=n_row, ncol=n_col )
        } else {
            x <- partable[ index_pp, "value"]
        }
        parm_list[[ pp ]] <- x
    }
    #--- output
    res <- list(partable=partable, parm_list=parm_list)
    return(res)
}
