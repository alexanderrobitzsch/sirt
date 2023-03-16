## File Name: dmlavaan_se_bootstrap.R
## File Version: 0.135

dmlavaan_se_bootstrap <- function(mod1, mod2, partable, R)
{
    parnames1 <- mod1$parnames
    NP1 <- mod1$NP
    parnames2 <- mod2$parnames
    NP2 <- mod2$NP
    est_boot1 <- matrix(NA, nrow=R, ncol=NP1)
    colnames(est_boot1) <- parnames1
    est_boot2 <- matrix(NA, nrow=R, ncol=NP2)
    colnames(est_boot2) <- parnames2
    partable1 <- mod1$partable
    partable1$start <- partable1$est
    partable2 <- mod2$partable
    partable2$start <- partable2$est

    #- reestimate models
    N <- mod1$N
    data <- mod1$args$data
    args1a <- mod1$args
    args1a$model <- partable1
    args2a <- mod2$args
    args2a$model <- partable2
    boot_sample <- matrix(0, nrow=N, ncol=R)
    rr <- 1
    while (rr<=R){
        # bootstrap data
        ind <- sort(sample(1:N, N, replace=TRUE))
        t1 <- table(ind)
        boot_sample[ as.numeric( names(t1) ), rr ] <- as.numeric(t1)
        # half sampling
        # ind <- sample(1:N, N/2)
        data_boot <- data[ind,]
        args1a$data <- data_boot
        args2a$data <- data_boot
        # estimate models
        mod1a <- try( do.call(what=mod1$fun, args=args1a), silent=TRUE)
        mod2a <- try( do.call(what=mod2$fun, args=args2a), silent=TRUE)
        accept <- ( ! inherits(mod1a,'try-error')) & ( ! inherits(mod1a,'try-error'))
        if (accept){
            rr <- rr + 1
            est_boot1[rr-1,1:NP1] <- coef(mod1a)
            est_boot2[rr-1,1:NP2] <- coef(mod2a)
        }
    }
    est_boot1 <- dmlavaan_remove_duplicated_columns(x=est_boot1)
    est_boot2 <- dmlavaan_remove_duplicated_columns(x=est_boot2)

    #-- now compute standard error for inference
    # pars_sel <- which( ! is.na( partable$diff ) )
    pars_sel <- 1:nrow(partable)
    pp <- pars_sel[1]
    for (pp in pars_sel){
        partable_pp <- partable[pp,]
        parname_pp <- partable_pp$parname
        ind_pp1 <- which( colnames(est_boot1)==parname_pp)[1]
        ind_pp2 <- which( colnames(est_boot2)==parname_pp)[1]
        if (length(ind_pp1)>0){
            partable$se1[pp] <- stats::sd(est_boot1[,ind_pp1])
        }
        if (length(ind_pp2)>0){
            partable$se2[pp] <- stats::sd(est_boot2[,ind_pp2])
        }
        if (length(ind_pp1)*length(ind_pp2)>0){
            z <- est_boot1[,ind_pp1]-est_boot2[,ind_pp2]
            partable$se_diff[pp] <- stats::sd(z)
        }
    }
    partable$t <- partable$diff / partable$se_diff
    partable$p <- 2*stats::pnorm(-abs(partable$t))

    #-- create bootstrap output objects
    res <- dmlavaan_se_bootstrap_create_est_boot(est_boot1=est_boot1,
                    est_boot2=est_boot2)
    est_boot <- res$est_boot
    V <- res$V

    #--- output
    res <- list(partable=partable, V=V, est_boot=est_boot,
                    boot_sample=boot_sample)
    return(res)
}
