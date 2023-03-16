## File Name: mgsem_proc_model_single_group.R
## File Version: 0.07

mgsem_proc_model_single_group <- function(model)
{
    H <- length(model)
    if (H==1){
        entries <- c("est","index")
        group1 <- list()
        group0 <- model[[1]]
        for (ee in entries){
            v1 <- list()
            for (vv in names(group0[[ee]]) ){
                v1[[vv]] <- 0*group0[[ee]][[vv]]
            } # end vv
            group1[[ee]] <- v1
        } # end ee
        res <- list( group0=group0, group1=group1)
    } else {
        res <- model
    }
    #-- output
    return(res)
}
