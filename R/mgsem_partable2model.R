## File Name: mgsem_partable2model.R
## File Version: 0.12


mgsem_partable2model <- function(partable, model, index=FALSE)
{
    ND <- nrow(partable)
    entries <- c("est")
    if (index){
        entries <- c("est","index")
    }
    for (entry in entries){
        for (dd in 1:ND){
            hh <- partable[dd,"group"]+1
            type <- paste(partable[dd,"type"])
            mat <- model[[hh]][[entry]][[type]]
            mat[ partable[dd,"i1"], partable[dd,"i2"] ] <- partable[dd,entry]
            model[[hh]][[entry]][[ type ]] <- mat
        }
    }
    return(model)
}
