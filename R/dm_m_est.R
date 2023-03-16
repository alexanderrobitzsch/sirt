## File Name: dm_m_est.R
## File Version: 0.106
## File Last Change: 2023-03-11

dm_m_est <- function(mod1, mod2)
{
    #*** create joint parameter table
    res <- dmlavaan_joint_parameterTable(mod1=mod1, mod2=mod2,
                    label_parnames='parnames')
    partable <- res$partable
    parnames <- res$parnames
    NP <- res$NP

    #*** sandwich estimate
    res <- dmlavaan_se_sandwich(mod1=mod1, mod2=mod2, partable=partable,
                    label_parnames='parnames', label_NPU='NP', label_B='B',
                    is_dmlavaan=TRUE)
    partable <- res$partable
    V <- res$partable

    #--- output
    res <- list(partable=partable, V=V, parnames=parnames, NP=NP)
    return(res)
}
