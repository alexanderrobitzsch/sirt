## File Name: dmlavaan_est_model.R
## File Version: 0.293

dmlavaan_est_model <- function(fun, args, h=1e-5, use_observed_info_lavaan=FALSE,
            method="sandwich")
{
    requireNamespace('lavaan')
    requireNamespace('MASS')
    remove_duplicated <- TRUE

    #- reorder groups
    if ('group' %in% names(args) ){
        data <- args$data
        group <- args$group
        data <- data[ order(data[,group]), ]
        args$data <- data
    }

    #- estimate model
    mod <- do.call(what=fun, args=args)
    data <- args$data
    N <- nrow(data)
    coef1 <- coef(mod)
    parnames <- names(coef1)
    NP <- length(coef1)
    #- standard errors
    vcov1 <- vcov(mod)

    #- parameter table
    res <- dmlavaan_est_model_parameterTable(mod=mod, parnames=parnames,
                    coef1=coef1, vcov1=vcov1)
    partable <- res$partable
    parnames <- res$parnames
    parnames0 <- res$parnames0
    coef1 <- res$coef1
    vcov1 <- res$vcov1
    NPU <- res$NPU
    se1 <- sqrt(diag(vcov1))

    #- score function
    if (method=='sandwich'){
        scores <- lavaan::lavScores(object=mod, remove.duplicated=remove_duplicated)
        colnames(scores) <- parnames0
        #--- compute derivative of score function
        use_observed_info_score_deriv <- TRUE
        if (use_observed_info_score_deriv){
            B <- dmlavaan_est_model_bread_matrix_score_derivatives(fun=fun, args=args,
                        partable=partable, scores=scores, h=h, symmetric=FALSE,
                        remove_duplicated=remove_duplicated, parnames=parnames0)
        }

        #- observed information matrix
        use_observed_info_lavaan <- FALSE
        if (use_observed_info_lavaan){
            BI <- lavaan::lavInspect( object=mod, what='information.observed')
            BI <- as.matrix(BI)
            sel <- partable[ partable$unique==1, 'free' ]
            B <- BI[sel,sel]
        }

        A <- crossprod(scores) / N^2
        res <- dmlavaan_sandwich_formula(A=A, B=B, parnames=parnames0)
        V <- res$V
        V_sw <- V
        se_sw <- res$se_sw

        #-- include results in partable
        partable <- dmlavaan_est_model_include_partable(partable=partable, se_sw=se_sw)
    } else {
        partable$se_sw <- partable$se
        A <- NULL
        B <- NULL
        scores <- NULL
        se_sw <- NULL
        V_sw <- NULL
    }

    #-- output
    res <- list(mod=mod, coef=coef1, vcov=vcov1, parnames=parnames, NP=NP,
                    se=se1, se_sw=se_sw, scores=scores, partable=partable,
                    V_sw=V_sw, A=A, B=B, N=N, NPU=NPU, data=data, fun=fun, args=args,
                    parnames0=parnames0)
    return(res)
}
