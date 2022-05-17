## File Name: meas_inv_cfa_modify_partable.R
## File Version: 0.136

meas_inv_cfa_modify_partable <- function(partable, mimod, critval)
{
    pars_mi <- meas_inv_compute_lavaan_parnames(object=mimod)
    pars_pt <- meas_inv_compute_lavaan_parnames(object=partable)
    free_parameter <- ( mimod$mi[1] > critval )
    if (free_parameter){
        ind <- which( pars_pt==pars_mi[1] & partable$group==mimod$group[1] )
        joint_label <- partable[ind,"label"]
        plabel <- partable[ind,"plabel"]
        is_joint_label <- joint_label==plabel
        partable[ind,"label"] <- ""
        if (is_joint_label){
            ind0 <- setdiff(which(pars_pt==pars_mi[1]),ind)
            partable0 <- partable[ind0, ]
            partable0 <- partable0[ partable0$label !="", ]
            plabel_new <- partable0$plabel[1]
            partable0[,"label"] <- plabel_new
            partable[ partable0$id, "label" ] <- partable0$label
            ind2 <- which( partable$op=="==" & partable$lhs==plabel)
            partable[ind2,"lhs"] <- plabel_new

            ind2 <- which( partable$op=="==" & (partable$lhs==partable$rhs) )
            if (length(ind2)>0){
                partable <- partable[-ind2, ]
            }

        }

        ind2 <- which( partable$op=="==" & (partable$lhs==plabel | partable$rhs==plabel))
        if (!is_joint_label){
            if (length(ind2)>0){
                partable <- partable[ - ind2, ]
            }
        }
    }
    #--- output
    res <- list(partable=partable, free_parameter=free_parameter, par_changed=pars_mi[1])
    return(res)
}


