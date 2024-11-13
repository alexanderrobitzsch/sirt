## File Name: linking.haberman.lq.R
## File Version: 0.255

linking.haberman.lq <- function(itempars, pow=2, eps=1e-3, a_log=TRUE,
    use_nu=FALSE, est_pow=FALSE, lower_pow=.1, upper_pow=3, method="joint",
    le=FALSE, vcov_list=NULL)
{
    CALL <- match.call()
    s1 <- Sys.time()

    use_pw <- method %in% c('pw1','pw2')
    if ( (!use_pw) | (pow!=2) | (!a_log) ){ le <- FALSE }

    #--- process item parameters
    res <- linking_proc_itempars(itempars=itempars)
    itempars <- res$itempars
    NS <- res$NS
    NI <- res$NI
    items <- res$items
    studies <- res$studies
    wgtM <- res$wgtM
    aM <- res$aM
    bM <- res$bM
    est_pars <- res$est_pars
    weights_exist <- res$weights_exist
    a.orig <- aM
    b.orig <- bM
    itempars1 <- itempars
    items <- sirt_sup(x=itempars1[,2])
    studies <- sirt_sup(x=itempars1[,1])
    G <- length(studies)
    I <- length(items)
    ind_items <- match(itempars1[,2], items)
    ind_studies <- match(itempars1[,1], studies)
    IP <- nrow(itempars1)
    X0 <- matrix(0, nrow=IP, ncol=G-1+I)
    colnames(X0) <- c(studies[-1], items)
    w <- itempars1[,5]

    res_optim <- list()
    eps_vec <- 10^seq(0,-10, by=-.5)
    item <- data.frame(item=items)

    #**** estimation of A
    y <- itempars1[,3]
    if (a_log){
        y <- log(y)
    }
    X <- X0
    for (gg in 2L:G){
        X[ ind_studies==gg, gg-1] <- 1
    }
    for (ii in 1L:I){
        X[ ind_items==ii, ii+G-1] <- 1
    }

    if ( use_pw ){
        res <- linking_haberman_lq_pw_create_design(y=y, ind_studies=ind_studies,
                    ind_items=ind_items, method=method)
        y <- res$y
        X <- res$X
        w <- res$w
        des_pw_slopes <- res
    }

    #- fit
    res_optim$slopes <- mod0 <- lq_fit(y=y, X=X, w=w, pow=pow, eps=eps,
                eps_vec=eps_vec, est_pow=est_pow, lower_pow=lower_pow,
                upper_pow=upper_pow)
    coef0 <- mod0$coefficients
    pow_slopes <- mod0$pow
    ind_groups <- 1L:(G-1)
    coef0_A <- coef0[ind_groups]
    if (!use_pw){
        a_joint <- coef0[-c(ind_groups)]
    } else {
        a_joint <- NULL
    }
    if (!use_pw){
        ar <- y - X %*% coef0
    } else {
        ar <- NULL
    }

    if (a_log){
        coef0_A <- exp(c(0,coef0_A))
        if (!use_pw){
            a_joint <- exp(a_joint)
            ar <- a_joint[ind_items]*( exp(ar) - 1 )
        }
    } else {
        coef0_A <- 1+c(0,coef0_A)
    }
    resid <- itempars1
    if (!use_pw){
        resid$adif <- ar
    }

    At <- coef0_A
    transf.personpars <- data.frame(study=studies, A_theta=coef0_A, se_A_theta=NA)
    item$a <- a_joint

    #**** estimation of B
    y <- itempars1[,4] * coef0_A[ ind_studies ]
    if (use_nu){
        y <- -itempars1[,4]*itempars1[,3]
        for (gg in 2L:G){
            ind_gg <- which(ind_studies==gg)
            X[ ind_gg, gg-1] <- itempars1[ ind_gg,3]/coef0_A[gg]
        }
    }
    if ( use_pw ){
        res <- linking_haberman_lq_pw_create_design(y=y, ind_studies=ind_studies,
                    ind_items=ind_items, method=method)
        y <- res$y
        X <- res$X
        w <- res$w
    }
    #- fit
    res_optim$intercepts <- mod0 <- lq_fit(y=y, X=X, w=w, pow=pow, eps=eps,
                eps_vec=eps_vec, est_pow=est_pow, lower_pow=lower_pow,
                upper_pow=upper_pow)
    coef0 <- mod0$coefficients

    pow_intercepts <- mod0$pow
    coef0_B <- -coef0[ind_groups]

    if (!use_pw){
        b_joint <- coef0[-c(ind_groups)]
    }
    if (use_nu){
        coef0_B <- -coef0_B
        if (!use_pw){
            b_joint <- -b_joint
        }
    }
    if (!use_pw){
        resid$b_resid <- y - X %*% coef0
    }
    Bt <- coef0_B <- c(0, coef0_B)
    transf.personpars$B_theta <- coef0_B
    transf.personpars$se_B_theta <- NA
    rownames(transf.personpars) <- NULL
    if (!use_pw){
        item$b <- b_joint
    }

    res_vcov <- NULL
    if (le){
        res_vcov <- linking_haberman_lq_pw_le(des=des_pw_slopes, res_optim=res_optim,
                        vcov_list=vcov_list)
    }  # end if le==TRUE

    #- include joint item parameters
    if (!use_pw){
        resid$a_joint <- item[ind_items, 'a']
        resid$b_joint <- item[ind_items, 'b']
    }

    # transformation for item parameters
    transf.itempars <- data.frame( study=studies, At=1/At, se_At=NA,
                            Bt=-Bt, se_Bt=NA )

    converged <- res_optim$intercepts$converged & res_optim$slopes$converged
    s2 <- Sys.time()
    time <- list(s1=s1, s2=s2)

    #-- output list
    description <- 'Linking based on L_q distance according to Haberman (2009)'
    res <- list( transf.personpars=transf.personpars, transf.itempars=transf.itempars,
                    pow=pow, eps=eps, item=item, resid=resid, description=description,
                    converged=converged, a_log=a_log, use_nu=use_nu, est_pow=est_pow,
                    pow_slopes=pow_slopes, pow_intercepts=pow_intercepts,
                    vcov=res_vcov, method=method, CALL=CALL, time=time)
    class(res) <- 'linking.haberman.lq'
    return(res)

}
