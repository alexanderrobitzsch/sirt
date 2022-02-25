## File Name: mgsem_proc_model_partable_define_index.R
## File Version: 0.11


mgsem_proc_model_partable_define_index <- function(partable)
{
    dfr <- partable
    ind1 <- which(dfr$index==1)
    N1 <- length(ind1)
    M1 <- max(dfr$index)+10
    if (length(ind1)>0){
        dfr$index[ind1] <- M1+1:N1
    }
    dfr$index <- match(dfr$index, unique(dfr$index))
    dfr$unique <- 1 * (! duplicated(dfr$index) )
    #** find duplicated parameters
    ND <- nrow(dfr)
    for (dd in 1:ND){
        if (dfr[dd,"unique"]==0){
            i1 <- which( dfr[,"index"]==dfr[dd,"index"])
            i2 <- which( dfr[,"unique"]==1)
            i3 <- intersect(i1,i2)
            dfr[dd,"name"] <- dfr[i3,"name"]
        }
        #- recycled parameters
        group_dd <- dfr[dd,"group"]
        if (group_dd>0){
            name_dd <- dfr[dd,"name"]
            name0_dd <- gsub(paste0("_G",group_dd),"_G0", name_dd)
            ind <- match(name0_dd, dfr$name)
            if (length(ind)==1){
                dfr[dd,"recycle"] <- ind
            }
        }
    }
    dfr$recycle <- ifelse(is.na(dfr$recycle), 0, dfr$recycle )

    dfr <- data.frame("rid"=1:ND, dfr)
    NP <- max(dfr$index)
    #--- output
    res <- list(partable=dfr, ND=ND, NP=NP)
    return(res)
}
