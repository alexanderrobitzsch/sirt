## File Name: mgsem_proc_model_update_penalties_matrix.R
## File Version: 0.073

mgsem_proc_model_update_penalties_matrix <- function(partable, entries, model)
{

    ND <- nrow(partable)
    for (entry in entries){
        for (dd in 1L:ND){
            if (partable[dd,'unique']==1){
                group <- partable$group[dd]+1
                type_dd <- paste(partable$type[dd])
                mat_gg_ee <- model[[group]][[entry]]
                if (!is.null(mat_gg_ee)){
                    mat_gg_ee <- mat_gg_ee[[type_dd]]
                    if (!is.null(mat_gg_ee)){
                        i1 <- partable$i1[dd]
                        i2 <- partable$i2[dd]
                        mat_gg_ee[ i1, i2] <- partable[dd,entry]
                        model[[group]][[entry]][[type_dd]] <- mat_gg_ee
                    }
                }
            }  # end is entry to be changed
        }  # end dd
    }  # end entry
    return(model)
}
