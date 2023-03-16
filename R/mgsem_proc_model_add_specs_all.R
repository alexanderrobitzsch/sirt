## File Name: mgsem_proc_model_add_specs_all.R
## File Version: 0.096
## File Last Change: 2022-05-16


mgsem_proc_model_add_specs_all <- function(model, entries, type, ii, jj, dfr1,
        group, N_group, names_prior_list=NULL, pen_type="scad")
{
    # N1 <- N_group[group]
    # if (group==0){
    #     N1 <- sum(N_group)
    # }
    N1 <- 1
    NE <- length(entries)
    for (entry in entries){
        if (entry=="lower"){ default <- -Inf }
        if (entry=="upper"){ default <- Inf }
        if (entry=="prior"){ default <- "none" }
        if (entry=="pen_l2"){ default <- 0 }
        if (entry=="pen_lp"){ default <- 0 }
        if (entry=="pen_difflp"){ default <- 0 }

        val <- mgsem_proc_model_add_specs(model=model, entry=entry, type=type,
                    ii=ii, jj=jj, default=default)
        if (entry=="prior"){ val <- paste(val) }
        if (entry %in% c("pen_l2","pen_lp","pen_difflp")){
            val <- N1*val
            if ( (entry=="pen_difflp") & (group==0)){
                val <- 0
            }
            #- no penalty
            if (pen_type=="none"){
                val <- 0
            }
        }

        #- correctness checks
        if (entry=="prior"){
            if ( ( val!="none") & ( ! val %in% names_prior_list ) ){
                v1 <- paste0( "Specified prior in ", type, "[",
                        ii, ",", jj, "] in group ", group, " not in 'prior_list'!\n")
                stop(v1)
            }
        }
        dfr1[1,entry] <- val
    }
    return(dfr1)
}
