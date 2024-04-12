## File Name: xxirt_createItemList.R
## File Version: 0.141


#--- create item list
xxirt_createItemList <- function( customItems, itemtype,items, partable )
{
    #*** list with customized items
    I <- length(items)
    item_list <- as.list( 1L:I )
    names(item_list) <- items
    CI <- length(customItems)
    for (ii in 1L:I){
        itemtype_ii <- itemtype[ii]
        item_list[[ii]] <- NA
        partable_ii <- partable[ partable$itemnr==ii, ]
        for (vv in 1L:CI){
            if ( paste(customItems[[vv]]$name)==itemtype_ii ){
                item_list[[ii]] <- customItems[[vv]]
            }  # end if
        } # end vv
        item_list[[ii]]$est <- partable_ii$est
        item_list[[ii]]$par <- partable_ii$value
        item_list[[ii]]$prior <- partable_ii$prior
        item_list[[ii]]$prior_par1 <- partable_ii$prior_par1
        item_list[[ii]]$prior_par2 <- partable_ii$prior_par2
    }  # end ii
    return(item_list)
}
