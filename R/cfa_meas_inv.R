## File Name: cfa_meas_inv.R
## File Version: 0.158


cfa_meas_inv <- function(dat, group, weights=NULL, alpha=0.01, verbose=FALSE,
                    op=c("~1","=~") )
{
    requireNamespace("lavaan")
    #--- define data input
    items <- colnames(dat)
    if (is.null(weights)){
        weights <- rep(1, nrow(dat))
    }
    dat <- data.frame(group=group, weight=weights, dat)
    args <- list(data=dat, group="group")
    args$sampling.weights <- "weight"

    #--- define model under scalar invariance
    lavmodel <- invariance_alignment_cfa_config_estimate_define_lavaan_model(
                        items_gg=items, label_F="F", model="2PM")
    args$model <- lavmodel
    args$meanstructure <- TRUE
    args$std.lv <- TRUE
    args$group.equal <- c("loadings","intercepts")

    mod1 <- lavaan::cfa(data=dat, model=lavmodel, std.lv=TRUE, meanstructure=TRUE,
                    group.equal=c("loadings","intercepts"), group="group",
                    sampling.weights="weight")
    partable1 <- lavaan::parameterTable(object=mod1)
    mimod1 <- mi_inv_lavaan_modification_indices(mod=mod1, op=op)
    partable0 <- partable1
    partable <- partable1
    mod0 <- mod1

    pars_mi <- meas_inv_cfa_proc_partable(partable=partable1, items=items)

    #--- change model by subsequently freeing parameters
    nfp <- 0
    critval <- stats::qchisq(1-alpha, df=1)
    free <- TRUE
    while(free){

        res <- meas_inv_cfa_modify_partable(partable=partable1, mimod=mimod1,
                            critval=critval)
        free <- res$free_parameter

        if (free){
            nfp <- nfp+1
            partable2 <- res$partable
            mod2 <- lavaan::cfa(data=dat, model=partable2, group="group",
                        sampling.weights="weight")
            partable1 <- lavaan::parameterTable(object=mod2)
            mimod1 <- mi_inv_lavaan_modification_indices(mod=mod2, op=op)
            if (verbose){
                cat(paste0("freed ", nfp, " parameters\n" ) )
                utils::flush.console()
            }

        }
        mod1 <- mod2
    }

    pars_pi <- meas_inv_cfa_proc_partable(partable=partable1, items=items)

    #--- output
    res <- list(pars_mi=pars_mi, pars_pi=pars_pi, alpha=alpha, critval=critval,
                    nfp=nfp, partable=partable1, dat=dat, items=items,
                    mod_mi=mod0, mod_pi=mod1)
    return(res)
}
