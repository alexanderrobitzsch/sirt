## File Name: rasch_mml2_mstep_one_step.R
## File Version: 1.071


rasch_mml2_mstep_one_step <- function(args0, prob_fun, entry, n.ik,
        diffindex, max.increment, numdiff.parm)
{
    h <- numdiff.parm
    val0 <- args0[[entry]]
    pjk <- do.call(what=prob_fun, args=args0)
    args1 <- rasch_mml2_modify_list_element( x=args0, entry=entry, value=val0+h )
    pjk1 <- do.call(what=prob_fun, args=args1)
    args1 <- rasch_mml2_modify_list_element( x=args0, entry=entry, value=val0-h )
    pjk2 <- do.call(what=prob_fun, args=args1)

    # numerical differentiation
    res <- rasch_mml2_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik,
                    diffindex=diffindex, max.increment=max.increment,
                    numdiff.parm=numdiff.parm )
    # update
    args0[[entry]] <- args0[[entry]] + res$increment
    res$args0 <- args0
    #- output
    return(res)
}
