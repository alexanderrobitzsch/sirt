## File Name: dmlavaan_est_model_parameterTable.R
## File Version: 0.072

dmlavaan_est_model_parameterTable <- function(mod, parnames, coef1, vcov1)
{
    partable <- lavaan::parameterTable(object=mod)
    partable$parname <- ''
    ind <- which( partable$free > 0)
    partable$parname[ind] <- parnames
    x1 <- setdiff( unique(partable$parname), '' )
    partable$pnid <- match( partable$parname, x1 )
    partable$pnid[ is.na(partable$pnid) ] <- 0
    partable$unique <- 1-duplicated(partable$parname)
    partable$unique[ partable$free==0 ] <- 0
    NPU <- max(partable$pnid)
    G <- max(partable$group)

    #* search for lavaan-defined parameter names
    pn1 <- paste(partable$parname)
    ld <- which( substring(pn1,1,2)=='.p')
    partable[ld, 'parname'] <- paste0( partable$lhs, partable$op, partable$rhs,
                                        '.g', partable$group )[ld]
    for (pp in 1:NPU){
        ind_pp <- which(partable$pnid==pp)
        lab1 <- partable[ind_pp[1],'parname']
        # if (length(ind_pp)>1){
        #     lab1 <- gsub('.g1', '', lab1, fixed=TRUE)
        # }
        partable[ind_pp,'parname'] <- lab1
    }

    parnames0 <- partable[partable$unique==1, 'parname']
    parnames <- partable[partable$free>0, 'parname']
    names(coef1) <- parnames
    colnames(vcov1) <- rownames(vcov1) <- parnames

    #--- output
    res <- list(partable=partable, parnames=parnames, G=G, NPU=NPU,
                    parnames0=parnames0, coef1=coef1, vcov=vcov)
    return(res)
}
