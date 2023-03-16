## File Name: linking.haberman.lq.R
## File Version: 0.198
## File Last Change: 2021-09-22

linking.haberman.lq <- function(itempars, pow=2, eps=1e-3, a_log=TRUE,
    use_nu=FALSE, est_pow=FALSE, lower_pow=.1, upper_pow=3)
{
    CALL <- match.call()
    s1 <- Sys.time()

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
    for (gg in 2:G){
        X[ ind_studies==gg, gg-1] <- 1
    }
    for (ii in 1:I){
        X[ ind_items==ii, ii+G-1] <- 1
    }
    #- fit
    res_optim$slopes <- mod0 <- lq_fit(y=y, X=X, w=w, pow=pow, eps=eps,
                eps_vec=eps_vec, est_pow=est_pow, lower_pow=lower_pow,
                upper_pow=upper_pow)
    coef0 <- mod0$coefficients
    pow_slopes <- mod0$pow
    ind_groups <- 1:(G-1)
    coef0_A <- coef0[ind_groups]
    a_joint <- coef0[-c(ind_groups)]
    ar <- y - X %*% coef0

    if (a_log){
        coef0_A <- exp(c(0,coef0_A))
        a_joint <- exp(a_joint)
        ar <- a_joint[ind_items]*( exp(ar) - 1 )
    } else {
        coef0_A <- 1+c(0,coef0_A)
    }
    resid <- itempars1
    resid$adif <- ar

    At <- coef0_A
    transf.personpars <- data.frame(study=studies, A_theta=coef0_A, se_A_theta=NA)
    item$a <- a_joint

    #**** estimation of B
    y <- itempars1[,4] * coef0_A[ ind_studies ]
    if (use_nu){
        y <- -itempars1[,4]*itempars1[,3]
        for (gg in 2:G){
            ind_gg <- which(ind_studies==gg)
            X[ ind_gg, gg-1] <- itempars1[ ind_gg,3]/coef0_A[gg]
        }
    }
    #- fit
    res_optim$intercepts <- mod0 <- lq_fit(y=y, X=X, w=w, pow=pow, eps=eps,
                eps_vec=eps_vec, est_pow=est_pow, lower_pow=lower_pow,
                upper_pow=upper_pow)
    coef0 <- mod0$coefficients
    pow_intercepts <- mod0$pow
    coef0_B <- -coef0[ind_groups]
    b_joint <- coef0[-c(ind_groups)]
    if (use_nu){
        coef0_B <- -coef0_B
        b_joint <- -b_joint
    }
    resid$b_resid <- y - X %*% coef0
    Bt <- coef0_B <- c(0, coef0_B)
    transf.personpars$B_theta <- coef0_B
    transf.personpars$se_B_theta <- NA
    rownames(transf.personpars) <- NULL
    item$b <- b_joint

    #- include joint item parameters
    resid$a_joint <- item[ind_items, "a"]
    resid$b_joint <- item[ind_items, "b"]

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
                    CALL=CALL, time=time)
    class(res) <- "linking.haberman.lq"
    return(res)

}
