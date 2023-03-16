## File Name: mgsem_observed_information.R
## File Version: 0.13


mgsem_observed_information <- function(coef, opt_fun_args, technical, comp_se=TRUE)
{
    if (technical$estimator!="ML"){
        comp_se <- FALSE
    }
    res <- list(coef=coef, comp_se=comp_se, se=NA*coef,
                    vcov=NULL, info_loglike=NULL, info_loglike_pen=NULL)

    if (comp_se){
        requireNamespace("MASS")
        grad_fun <- 'mgsem_numerical_gradient'
        h <- technical$h

        #- observed information without penalty
        opt_fun_args1 <- opt_fun_args
        opt_fun_args1$use_penalty <- FALSE
        args_grad <- list(par=coef, FUN=mgsem_grad_fun, symmetrize=TRUE,
                    h=h, opt_fun_args=opt_fun_args1, output_all=FALSE)
        res$info_loglike <- do.call(what=grad_fun, args=args_grad)

        #- observed information with penalty
        args_grad$opt_fun_args <- opt_fun_args
        res$info_loglike_pen <- do.call(what=grad_fun, args=args_grad)

        #- standard errors
        res$vcov <- MASS::ginv(X=res$info_loglike)
        res$se <- sqrt(diag(res$vcov))

    }

    #-- output
    return(res)
}
