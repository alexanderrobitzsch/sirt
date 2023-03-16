## File Name: mgsem_modify_model.R
## File Version: 0.04
## File Last Change: 2022-04-19


mgsem_modify_model <- function(model, group, entry, type, value)
{
    for (gg in group){
        for (ee in entry){
            for (tt in type){
                model[[gg+1]][[ee]][[tt]] <- value+0*model[[gg+1]][[ee]][[tt]]
            }  # end tt
        }  # end ee
    }  # end gg
    #--- output
    return(model)
}
