## File Name: lsem_fitsem_joint_estimation_prepare_partable_include_group_label.R
## File Version: 0.10


lsem_fitsem_joint_estimation_prepare_partable_include_group_label <- function(
    partable, gg, label_list)
{
    hh <- "lhs"
    for (hh in c("lhs","rhs")){
        ind <- which( partable[,hh] %in% setdiff(label_list, "") )
        if (length(ind)>0){
            partable[ ind, hh] <- paste0(partable[ ind, hh], "g", gg)
        }
    }
    partable$plabel <- paste0(label_list,"g", gg)
    partable$plabel[ paste(label_list)=="" ] <- ""
    return(partable)
}
