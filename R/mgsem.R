## File Name: mgsem.R
## File Version: 0.497

mgsem <- function(suffstat, model, data=NULL, group=NULL, weights=NULL,
        estimator="ML", p_me=2, p_pen=1, pen_type="scad",
        diffpar_pen=NULL, a_scad=3.7, eps_approx=1e-3, comp_se=TRUE,
        se_delta_formula=FALSE,
        prior_list=NULL, hessian=TRUE,
        fixed_parms=FALSE, partable_start=NULL,
        num_approx=FALSE, technical=NULL, control=list() )
{
    #- pen_type: lasso, scad or none

    #*** preliminaries
    CALL <- match.call()
    s1 <- Sys.time()

    test <- FALSE
    # test <- TRUE

    #*** process data if suffstat not available
    if (missing(suffstat)){
        data_proc <- mgsem_proc_data(data=data, group=group, weights=weights)
        groups <- data_proc$groups
        suffstat <- data_proc$suffstat
        weights <- data_proc$weights
        is_data <- TRUE
        G <- length(groups)
    } else {
        groups <- names(suffstat)
        G <- length(groups)
        if (is.null(groups)){
            groups <- paste0('Group',1:G)
        }
        data_proc <- NULL
        is_data <- FALSE
    }
    technical$is_data <- is_data

    #*** compute covariance matrix of sufficient statistics
    suffstat_vcov <- mgsem_suffstat_covariance_matrix(suffstat=suffstat)

    #*** process technical defaults
    res <- mgsem_proc_technical(technical=technical, control=control, p_me=p_me,
                p_pen=p_pen, eps_approx=eps_approx, suffstat=suffstat,
                estimator=estimator, diffpar_pen=diffpar_pen)
    technical <- res$technical
    technical$estimator <- estimator
    control <- res$control
    eps_approx <- res$eps_approx
    p <- res$p

    #*** process sufficient statistics
    res <- mgsem_proc_suffstat(suffstat=suffstat)
    suffstat <- res$suffstat
    G <- res$G
    I <- res$I
    N <- res$N
    N_group <- res$N_group
    random_sd <- -9
    if (test){
        random_sd <- 1e-1
    }

    #*** process model specification
    res <- mgsem_proc_model(model=model, G=G, prior_list=prior_list,
                    technical=technical, N_group=N_group, random_sd=random_sd,
                    pen_type=pen_type, fixed_parms=fixed_parms,
                    partable_start=partable_start, diffpar_pen=diffpar_pen)
    model <- res$model
    partable <- res$partable
    NP <- res$NP
    coef <- res$coef
    ND <- res$ND
    D <- res$D
    is_B <- res$is_B
    technical <- res$technical
    types <- res$types
    difflp_info <- res$difflp_info
    loop_parms <- res$loop_parms

    if (estimator=='ML'){
        eval_fun <- 'mgsem_loglike_suffstat'
        grad_param_fun <- 'mgsem_loglike_suffstat_derivative_parameter'
        grad_suffstat_fun <- 'mgsem_loglike_suffstat_derivative'
    }
    if (estimator=='ME'){
        eval_fun <- 'mgsem_loss_function_suffstat'
        grad_param_fun <- 'mgsem_loss_function_suffstat_derivative_parameter'
        grad_suffstat_fun <- 'mgsem_loss_function_suffstat'
    }

    if (technical$use_deriv){
        grad_fun <- mgsem_grad_fun
    } else {
        grad_fun <- mgsem_grad_fun_numeric_approx
    }

    x <- coef
    opt_fun_args <- list(NP=NP, ND=ND, suffstat=suffstat, model=model,
                        partable=partable, G=G, I=I, D=D, is_B=is_B,
                        N=N, N_group=N_group,
                        estimator=estimator, eval_fun=eval_fun,
                        grad_param_fun=grad_param_fun,
                        grad_suffstat_fun=grad_suffstat_fun,
                        technical=technical, p_pen=p_pen, p_me=p_me,
                        p=p, eps_approx=eps_approx, num_approx=num_approx,
                        prior_list=prior_list, use_penalty=TRUE, types=types,
                        difflp_info=difflp_info, loop_parms=loop_parms,
                        pen_type=pen_type, a_scad=a_scad )

    #- test function
    res <- mgsem_test_fun(test=test, coef=coef, opt_fun_args=opt_fun_args)

    #**** estimation
    if (! fixed_parms){

        #- define lower and upper bounds
        ind <- which(partable$unique>0)
        lower <- partable[ ind, 'lower' ]
        upper <- partable[ ind, 'upper' ]

        #- use optimizer
        opt_res <- sirt_optimizer(optimizer=technical$optimizer, par=x, fn=mgsem_opt_fun,
                        grad=grad_fun, opt_fun_args=opt_fun_args,
                        method='L-BFGS-B', lower=lower, upper=upper,
                        hessian=hessian, control=control )
    } else {
        #**** no estimation
        ll <- mgsem_opt_fun(x=x, opt_fun_args=opt_fun_args)
        opt_res <- list(par=x, opt_value=ll)
    }

    #-- collect results
    coef <- opt_res$par
    opt_value <- opt_res$value
    partable <- mgsem_coef2partable(coef=coef, partable=partable)
    model <- mgsem_partable2model(partable=partable, model=model)

    #-- evaluate optimization function and gradient
    opt_fun_args$partable <- partable
    opt_fun_args$model <- model
    opt_fun_output <- mgsem_opt_fun(x=coef, opt_fun_args=opt_fun_args, output_all=TRUE)
    implied <- opt_fun_output$implied
    est_tot <- opt_fun_output$est_tot
    grad_fun_output <- mgsem_grad_fun(x=coef, opt_fun_args=opt_fun_args, output_all=TRUE)

    #-- vcov for estimator='ME'
    res <- mgsem_vcov_me(coef=coef, opt_fun_args=opt_fun_args,
                            suffstat_vcov=suffstat_vcov, comp_se=comp_se,
                            se_delta_formula=se_delta_formula)
    vcov <- res$vcov
    se <- res$se
    comp_se_me <- res$comp_se_me

    #-- residual statistics
    residuals <- mgsem_output_proc_residuals(implied=implied, suffstat=suffstat)

    #-- compute case-wise log-likelihood
    case_ll <- mgsem_output_proc_casewise_likelihood(data_proc=data_proc,
                    implied=implied, estimator=estimator)

    #-- computation of standard errors
    res <- mgsem_observed_information(coef=coef, opt_fun_args=opt_fun_args,
                    technical=technical, comp_se=comp_se, comp_se_me=comp_se_me)
    if (!comp_se_me){
        vcov <- res$vcov
        comp_se <- res$comp_se
        se <- res$se
    }
    info_loglike <- res$info_loglike
    info_loglike_pen <- res$info_loglike_pen

    #-- information criteria
    ic <- mgsem_ic(opt_fun_output=opt_fun_output, opt_fun_args=opt_fun_args,
                        partable=partable, technical=technical)

    #--- output
    s2 <- Sys.time()
    res <- list(coef=coef, vcov=vcov, se=se,
                    partable=partable, model=model, suffstat=suffstat,
                    opt_res=opt_res, opt_value=opt_value, implied=implied,
                    est_tot=est_tot, residuals=residuals, opt_fun_output=opt_fun_output,
                    ic=ic, info_loglike=info_loglike, info_loglike_pen=info_loglike_pen,
                    estimator=estimator, p_pen=p_pen, p_me=p_me,
                    pen_type=pen_type, eps_approx=eps_approx,
                    technical=technical, comp_se=comp_se, groups=groups, group=group,
                    data=data, data_proc=data_proc, case_ll=case_ll,
                    suffstat_vcov=suffstat_vcov, CALL=CALL, s1=s1, s2=s2)
    class(res) <- 'mgsem'
    return(res)
}
