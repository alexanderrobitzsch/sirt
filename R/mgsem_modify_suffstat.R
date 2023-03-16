## File Name: mgsem_modify_suffstat.R
## File Version: 0.01


mgsem_modify_suffstat <- function(model, group, entry, value)
{
    for (gg in group){
        for (ee in entry){
            model[[gg]][[ee]] <- value+0*model[[gg]][[ee]]
        }  # end ee
    }  # end gg
    #--- output
    return(model)
}
