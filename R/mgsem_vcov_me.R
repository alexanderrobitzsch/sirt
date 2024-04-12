## File Name: mgsem_vcov_me.R
## File Version: 0.171

mgsem_vcov_me <- function(coef, opt_fun_args, suffstat_vcov, comp_se,
            se_delta_formula=FALSE)
{
    estimator <- opt_fun_args$technical$estimator
    comp_se_me <- (estimator=='ME') & comp_se
    if ( (estimator=='ML') & comp_se & se_delta_formula ){
        comp_se_me <- TRUE
    }
    A <- NULL
    V <- NULL
    if (comp_se_me){
        NP <- length(coef)
        parnames <- names(coef)
        #-- derivative with respect to coef
        grad_der_coef <- matrix(0, nrow=NP, ncol=NP)
        rownames(grad_der_coef) <- colnames(grad_der_coef) <- parnames
        h <- opt_fun_args$technical$h
        for (pp in 1L:NP){
            coef1 <- mgsem_add_increment(x=coef, h=h, i1=pp)
            res1 <- mgsem_grad_fun(x=coef1, opt_fun_args=opt_fun_args, output_all=FALSE)
            coef2 <- mgsem_add_increment(x=coef, h=-h, i1=pp)
            res2 <- mgsem_grad_fun(x=coef2, opt_fun_args=opt_fun_args, output_all=FALSE)
            grad_der_coef[,pp] <- ( res1 - res2 )/(2*h)
        }
        grad_der_coef <- ( grad_der_coef + t(grad_der_coef) ) / 2

        #-- derivative with respect input parameters of sufficient statistics
        suffstat_pars <- suffstat_vcov$suffstat_pars
        SP <- nrow(suffstat_pars)
        V <- suffstat_vcov$suffstat_vcov
        rownames(V) <- suffstat_pars$label
        colnames(V) <- suffstat_pars$label
        grad_der_suffstat <- matrix(0, nrow=NP, ncol=SP)
        rownames(grad_der_suffstat) <- parnames
        colnames(grad_der_suffstat) <- suffstat_pars$label
        opt_fun_args1 <- opt_fun_args
        suffstat <- opt_fun_args$suffstat
        for (pp in 1L:SP){
            suffstat_pars_pp <- suffstat_pars[pp,]
            group_pp <- suffstat_pars_pp$group
            val_pp <- list(NA,2)
            for (oo in 1L:2){
                suffstat1 <- suffstat
                u <- h
                if (oo==2){
                    u <- -h
                }
                if (suffstat_pars_pp$type=='mu'){
                    entry1 <- mgsem_add_increment(x=suffstat1[[group_pp]]$M, h=u,
                                        i1=suffstat_pars_pp$index1)
                    suffstat1[[group_pp]]$M <- entry1
                } else {
                    entry1 <- mgsem_add_increment(x=suffstat1[[group_pp]]$S, h=u,
                                        i1=suffstat_pars_pp$index1,
                                        i2=suffstat_pars_pp$index2, symm=TRUE)
                    suffstat1[[group_pp]]$S <- entry1
                }
                opt_fun_args1$suffstat <- suffstat1
                val_pp[[oo]] <- mgsem_grad_fun(x=coef, opt_fun_args=opt_fun_args1,
                                    output_all=FALSE)
            }
            der_est <- (val_pp[[1]]-val_pp[[2]])/(2*h)

            grad_der_suffstat[,pp] <- der_est
        }

        #-- compute transformation matrix A
        W1 <- mgsem_ginv(X=grad_der_coef)
        A <- - ( W1 %*% grad_der_suffstat )
        vcov <- A %*% V %*% t(A)
        rownames(vcov) <- colnames(vcov) <- names(coef)
        se <- mgsem_sqrt_diag(x=vcov)
        names(se) <- names(coef)
    } else {
        vcov <- NULL
        se <- NULL
    }

    #--- output
    res <- list(vcov=vcov, se=se, A=A, V=V, comp_se_me=comp_se_me)
    return(res)
}
