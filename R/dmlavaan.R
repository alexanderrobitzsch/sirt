## File Name: dmlavaan.R
## File Version: 0.164
## File Last Change: 2023-03-10


dmlavaan <- function( fun1, args1, fun2, args2, method="sandwich", R=50)
{
    requireNamespace('lavaan')
    #-- estimate models
    mod1 <- dmlavaan_est_model(fun=fun1, args=args1, method=method)
    mod2 <- dmlavaan_est_model(fun=fun2, args=args2, method=method)

    #-- joint parameter table
    res <- dmlavaan_joint_parameterTable(mod1=mod1, mod2=mod2)
    partable <- res$partable
    NP <- res$NP
    parnames <- res$parnames

    #* NULL objects
    est_boot <- NULL
    boot_sample <- NULL
    V <- NULL

    #-- sandwich estimate
    if (method=='sandwich'){
        res <- dmlavaan_se_sandwich(mod1=mod1, mod2=mod2, partable=partable)
        partable <- res$partable
        V <- res$V
    }

    #-- bootstrap estimate
    if (method=='bootstrap'){
        res <- dmlavaan_se_bootstrap(mod1=mod1, mod2=mod2,
                            partable=partable, R=R)
        partable <- res$partable
        V <- res$V
        est_boot <- res$est_boot
        boot_sample <- res$boot_sample
    }

    #- create vector of coefficients
    coef <- dmlavaan_create_coef(partable=partable)

    #** define output
    partable1 <- mod1$partable
    partable2 <- mod2$partable

    #-- output
    res <- list(coef=coef, V=V, partable=partable, mod1=mod1, mod2=mod2, method=method,
                    parnames=parnames, NP=NP, partable1=partable1,
                    partable2=partable2, est_boot=est_boot,
                    boot_sample=boot_sample, R=R)
    return(res)
}
