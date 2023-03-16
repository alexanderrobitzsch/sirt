## File Name: dmlavaan_joint_parameterTable.R
## File Version: 0.06
## File Last Change: 2023-03-11

dmlavaan_joint_parameterTable <- function(mod1, mod2, label_parnames="parnames0")
{
    parnames1 <- mod1[[ label_parnames ]]
    parnames2 <- mod2[[ label_parnames ]]
    #-- create joint parameter table
    parnames <- union(parnames1, parnames2)
    NP <- length(parnames)
    partable <- data.frame(id=1:NP, parname=parnames)
    partable$in_mod1 <- match(partable$parname, parnames1)
    partable$in_mod2 <- match(partable$parname, parnames2)
    #- first model
    partable <- dmlavaan_joint_parameterTable_merge_table(partable=partable,
                        mod=mod1, model_index=1)
    #- second model
    partable <- dmlavaan_joint_parameterTable_merge_table(partable=partable,
                        mod=mod2, model_index=2)
    partable$diff <- partable$est1 - partable$est2
    partable$se_diff <- NA
    rownames(partable) <- NULL

    #--- output
    res <- list(partable=partable, parnames=parnames, NP=NP)
    return(res)
}
