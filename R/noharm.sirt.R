## File Name: noharm.sirt.R
## File Version: 0.925


########################################
# NOHARM implementation in R
########################################

#---------------
# NOHARM
# Normal-Ogive Harmonic Analysis Robust Method
#---------------
# The other function with linear and loglinear models of parameters
# is denoted as 'noharm.loglm.sirt'.
#---------------

noharm.sirt <- function(dat, pm=NULL, N=NULL, weights=NULL, Fval=NULL, Fpatt=NULL,
    Pval=NULL, Ppatt=NULL, Psival=NULL, Psipatt=NULL, dimensions=NULL,
    lower=0, upper=1, wgtm=NULL, pos.loading=FALSE, pos.variance=FALSE,
    pos.residcorr=FALSE, maxiter=1000, conv=1e-6, optimizer="nlminb",
    par_lower=NULL, reliability=FALSE, ... )
{
    v1 <- s1 <- Sys.time()
    CALL <- match.call()
    #*** data preprocessing
    if (missing(dat)){
        input_pm <- TRUE
        dat <- NULL
    } else {
        input_pm <- FALSE
    }
    res <- noharm_sirt_preproc( dat=dat, pm=pm, N=N, weights=weights, Fpatt=Fpatt,
                Fval=Fval, Ppatt=Ppatt, Pval=Pval, Psipatt=Psipatt, Psival=Psival,
                wgtm=wgtm, dimensions=dimensions, pos.loading=pos.loading,
                pos.variance=pos.variance, pos.residcorr=pos.residcorr,
                input_pm=input_pm, lower=lower, upper=upper)
    ss <- res$ss
    pm0 <- res$pm0
    sumwgtm <- res$sumwgtm
    model.type <- res$model.type
    modtype <- res$modtype
    N <- res$N
    I <- res$I
    dat0 <- res$dat0
    Fpatt <- res$Fpatt
    Fval <- res$Fval
    Ppatt <- res$Ppatt
    Pval <- res$Pval
    Psipatt <- res$Psipatt
    Psival <- res$Psival
    D <- res$D
    pm <- res$pm
    dat <- res$dat
    dat.resp <- res$dat.resp
    wgtm <- res$wgtm
    F_dimnames <- res$F_dimnames
    items <- res$items
    wgtm.default <- res$wgtm.default
    parm_table <- res$parm_table
    npar <- res$npar
    parm_index <- res$parm_index
    lower <- res$lower
    upper <- res$upper

    modesttype <- 1   # NOHARM estimation

    #----
    # compute betaj
    if (modesttype==1){
        betaj <- - stats::qnorm( ( diag(pm) - lower ) / ( upper - lower ) )
    }
    # compute bjk coefficients
    # include lower and upper asymptotes here

    b0 <- lower + (upper-lower) * stats::pnorm(-betaj)
    b1 <- (upper-lower) * sirt_dnorm(betaj)
    b2 <- betaj * b1 / sqrt(2)
    b3 <- ( betaj^2 -1 ) * b1 / sqrt(6)

    # create fixed cofficients
    b0_jk <- b0.jk <- noharm_sirt_outer_coefs(x=b0)
    b1_jk <- b1.jk <- noharm_sirt_outer_coefs(x=b1)
    b2_jk <- b2.jk <- noharm_sirt_outer_coefs(x=b2)
    b3_jk <- b3.jk <- noharm_sirt_outer_coefs(x=b3)

    par0 <- noharm_sirt_partable_extract_par(parm_table=parm_table)
    parm_table_free <- parm_table[ parm_table$fixed==0, ]

    args_optim <- list(parm_table=parm_table, parm_index=parm_index, I=I, D=D,
                    b0.jk=b0.jk, b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk,
                    pm=pm, wgtm=wgtm, use_rcpp=TRUE)

    optim_fn <- function(x){
        args_optim$x <- x
        val <- do.call(what=noharm_sirt_optim_function, args=args_optim)
        return(val)
    }
    optim_gr <- function(x){
        args_optim$x <- x
        grad <- do.call(what=noharm_sirt_optim_gradient, args=args_optim)
        return(grad)
    }

    v2 <- Sys.time()
    time <- list( time_pre=v2-v1)

    #* optimization
    if (is.null(par_lower)){
        par_lower <- noharm_sirt_partable_extract_par(parm_table=parm_table, col='lower')
    }
    res <- sirt_optimizer(optimizer=optimizer, par=par0, fn=optim_fn,
                grad=optim_gr, lower=par_lower, method='L-BFGS-B', hessian=FALSE, ...)
    res_opt <- res
    par0 <- res$par
    v3 <- Sys.time()
    time$time_opt <- v3-v2

    #* include in parameter table
    parm_table <- noharm_sirt_partable_include_par(par=par0, parm_table=parm_table)
    Fval <- Fmat <- noharm_sirt_create_parameter_matrices('F', parm_table=parm_table,
                                parm_index=parm_index)
    Pval <- Pmat <- noharm_sirt_create_parameter_matrices('P', parm_table=parm_table,
                                parm_index=parm_index)
    Psival <- Psimat <- noharm_sirt_create_parameter_matrices('Psi',
                                parm_table=parm_table, parm_index=parm_index)
    estpars <- list( estF=sum(Fpatt>0), estP=sum(Ppatt>0), estPsi=sum(Psipatt>0) )

    #--- include colnames in estimated matrices
    colnames(Fval) <- F_dimnames
    rownames(Fval) <- items
    rownames(Pval) <- colnames(Pval) <- F_dimnames
    rownames(Psival) <- colnames(Psival) <- items

    #* compute final constants
    res <- noharm_sirt_compute_final_constants( Fval=Fval, Pval=Pval, betaj=betaj,
                        modesttype=modesttype )
    f0 <- res$f0
    uqn <- res$uqn
    loadingsF <- res$loadingsF

    #* compute residuals
    res <- noharm_sirt_est_residuals( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I,
                D=D, b0.jk=b0.jk, b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk, wgtm=wgtm, pm=pm,
                Psival=Psival, Psipatt=Psipatt, modesttype=modesttype )
    residuals <- res$residuals
    rmsr <- res$rmsr
    tanaka <- res$tanaka

    #**** arrange output list
    res <- list( tanaka=tanaka, rmsr=rmsr, N.itempair=ss, pm=pm0, wgtm=wgtm, sumwgtm=sumwgtm,
                lower=lower, upper=upper, residuals=residuals, final.constants=f0,
                factor.cor=Pval, thresholds=-betaj, uniquenesses=uqn, loadings=loadingsF,
                loadings.theta=Fval, residcorr=Psival, model.type=model.type, modtype=modtype,
                Nobs=N, Nitems=I, Fpatt=Fpatt, Ppatt=Ppatt, Psipatt=Psipatt,
                dat=dat0, systime=Sys.time(), dimensions=D, display.fit=5,
                res_opt=res_opt, parm_table=parm_table, b0=b0, b1=b1, b2=b2, b3=b3,
                par_lower=par_lower)

    #*** number of estimated parameters
    Nestpars <- noharm_sirt_number_estimated_parameters( I=I, Fpatt=Fpatt,
                    Ppatt=Ppatt, Psipatt=Psipatt )
    res$Nestpars <- Nestpars

    #*** compute chi square statistics
    res <- noharm_sirt_compute_chi_square_statistics( res=res,
                residuals=residuals, pm=pm, N=N, I=I, sumwgtm=sumwgtm, modesttype=modesttype,
                Nestpars=Nestpars )

    #*** compute Green-Yang reliability
    v1 <- 1 * ( ( wgtm - diag(wgtm) ) > 0 )
    # ensure positive definiteness
    res$omega.rel <- NA
    L0 <- sqrt( diag(Pval) )
    N1 <- length(L0)
    L1inv <- L1 <- matrix( 0, nrow=N1, ncol=N1 )
    rownames(L1) <- rownames(L1inv) <- colnames(L1) <- colnames(L1inv) <- rownames(Pval)
    diag(L1) <- L0
    diag(L1inv) <- 1 / L0
    Fval2 <- Fval %*% L1
    Pval2 <- L1inv %*% Pval %*% L1inv
    dj <- sqrt( diag( Fval2 %*% Pval2 %*% t(Fval2) ) )
    ej <- sqrt( 1+dj^2 )
    Fval2 <- Fval2 / ej
    standardized.solution <- list( Fval=Fval2, Pval=Pval2 )
    if (reliability){
        res0 <- reliability.nonlinearSEM(facloadings=Fval2,
            thresh=res$thresholds, resid.cov=res$residcorr, cor.factors=Pval2 )
        res$omega.rel <- res0$omega.rel
    }
    if ( sum(v1) + I < I^2 ){
        res$omega.rel <- NA
    }

    #*** rotated solution for EFA
    if (model.type=='EFA'){
        res <- noharm_sirt_efa_rotated_solution( res=res, items=items,
                    F_dimnames=F_dimnames )
    }
    if (modesttype==2){
        res$tetracor <- pm
    }

    #***----- more output
    res$estpars <- estpars
    res$modesttype <- modesttype
    res$guesses <- res$lower
    res$wgtm.default <- wgtm.default
    res$standardized.solution <- standardized.solution
    s2 <- Sys.time()
    res$s1 <- s1
    res$s2 <- s2
    res$CALL <- CALL
    v4 <- Sys.time()
    time$time_post <- v4-v3
    res$time <- time
    class(res) <- 'noharm.sirt'
    return(res)
}
