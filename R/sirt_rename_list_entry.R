## File Name: sirt_rename_list_entry.R
## File Version: 0.09

sirt_rename_list_entry <- function(model, name1, name2)
{
    model_temp <- model
    G <- length(model)
    for (gg in 1:G){
        model_gg <- model_temp[[gg]]
        TP <- length(model_gg)
        for (tt in 1:TP){
            model_gg_tt <- model_gg[[tt]]
            ind <- which( names(model_gg_tt)==name1)
            if (length(ind)>0){
                names(model_gg_tt)[ind] <- name2
            }
            model_gg[[tt]] <- model_gg_tt
        }
        model_temp[[gg]] <- model_gg
    }
    return(model_temp)
}
