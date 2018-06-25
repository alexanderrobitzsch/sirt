## File Name: rm_sdt_create_partable.R
## File Version: 0.45


rm_sdt_create_partable <- function( item.index, rater.index,
    est.c.rater, est.d.rater, tau.item, c.rater, diffindex,
    tau.prior, a.prior, d.prior, c.prior, est.a.item )
{
    I <- nrow(tau.item)
    ND <- nrow(c.rater)
    K <- ncol(c.rater)

    #--- tau.item
    v1 <- 1:(K*I) + 0*tau.item
    partable <- rm_sdt_create_parm_index_modify_elements(x=v1, start_index=1, type="tau")

    #--- a.item
    v1 <- seq_len(I)
    if (! est.a.item){
        v1 <- -9 + 0*v1
    }
    start <- max(partable$parindex) + 1
    v1 <- rm_sdt_create_parm_index_modify_elements(x=v1, start_index=start, type="a" )
    partable <- rbind( partable, v1)
    partable_item <- partable

    #--- c.rater
    v1 <- 0*c.rater
    g1 <- rm_sdt_create_parm_index_rater( est.rater=est.c.rater, ND=ND,
                item.index=item.index, rater.index=rater.index )
    M <- 0
    for (kk in 1:K){
        v1[,kk] <- g1 + M
        M <- max(v1[,kk])
    }
    start <- 1
    partable <- rm_sdt_create_parm_index_modify_elements(x=v1, start_index=start, type="c")

    #--- d.rater
    v1 <- rm_sdt_create_parm_index_rater( est.rater=est.d.rater, ND=ND,
                item.index=item.index, rater.index=rater.index )
    start <- max(partable$parindex) + 1
    v1 <- rm_sdt_create_parm_index_modify_elements(x=v1, start_index=start, type="d")
    partable <- rbind( partable, v1)
    partable_rater <- partable

    partable_item <- rm_sdt_create_partable_include_index(partable=partable_item)
    partable_rater <- rm_sdt_create_partable_include_index(partable=partable_rater)

    #--- index lists
    par_index <- list()
    par_index$tau.item <- which(partable_item$type=="tau")
    par_index$a.item <- which(partable_item$type=="a")
    par_index$c.rater <- which(partable_rater$type=="c")
    par_index$d.rater <- which(partable_rater$type=="d")

    #--- parameter groups for differentiation
    partable_item <- rm_sdt_create_partable_define_pargroups(partable=partable_item,
                            pg1="tau", pg2="a")
    partable_rater <- rm_sdt_create_partable_define_pargroups(partable=partable_rater,
                            pg1="c", pg2="d")

    #--- diffindex pargroups
    pargroup_item <- rm_sdt_create_partable_pargroup_indices( partable=partable_item,
                    item.index=item.index, diffindex=diffindex )
    pargroup_rater <- rm_sdt_create_partable_pargroup_indices( partable=partable_rater,
                    item.index=item.index, diffindex=diffindex )

    #--- prior distributions
    partable$prior_M <- NA
    partable$prior_SD <- NA
    partable_item <- rm_sdt_create_partable_include_priors(partable=partable_item,
                                type="tau", obj=tau.prior)
    partable_item <- rm_sdt_create_partable_include_priors(partable=partable_item,
                                type="a", obj=a.prior)
    partable_rater <- rm_sdt_create_partable_include_priors(partable=partable_rater,
                                type="c", obj=c.prior)
    partable_rater <- rm_sdt_create_partable_include_priors(partable=partable_rater,
                                type="d", obj=d.prior)
    #--- output
    res <- list( partable_item=partable_item, partable_rater=partable_rater,
                par_index=par_index, pargroup_item=pargroup_item,
                pargroup_rater=pargroup_rater )
    return(res)
}
