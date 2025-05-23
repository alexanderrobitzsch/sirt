## File Name: invariance.alignment.R
## File Version: 4.019


invariance.alignment <- function( lambda, nu, wgt=NULL,
    align.scale=c(1,1), align.pow=c(.5,.5), eps=1e-3,
    psi0.init=NULL, alpha0.init=NULL, center=FALSE, optimizer="optim",
    fixed=NULL, meth=1, vcov=NULL, eps_grid=seq(0,-10, by=-.5),
    num_deriv=FALSE, ... )
{
    CALL <- match.call()
    s1 <- Sys.time()
    type <- 'AM'
    align.pow0 <- align.pow
    align.pow <- align.pow / 2
    overparam <- FALSE
    le <- FALSE

    #-- labels for groups and items
    lambda0 <- lambda <- invariance_alignment_proc_labels(x=lambda)
    nu0 <- nu <- invariance_alignment_proc_labels(x=nu)

    #-- weights
    G <- nrow(lambda)   # number of groups
    I <- ncol(lambda)   # number of items
    if ( is.null(wgt) ){
        wgt <- 1+0*nu
    }

    #- reparametrization of meth argument
    meth0 <- meth
    if (meth0==4){ meth <- 0}     # Mplus FREE
    if (meth0==3){ meth <- 0.5}      # Mplus FIXED

    #- choose fixed value
    fixed <- invariance_alignment_choose_fixed(fixed=fixed, G=G, Gmax=999)
    reparam <- ! fixed
    if (meth%in%c(0,0.5)){
        constraint <- 'prod'
        reparam <- TRUE
        num_deriv <- TRUE
    }
    if (overparam){
        num_deriv <- TRUE
    }

    W1 <- dim(wgt)
    wgtM <- matrix( colSums(wgt,na.rm=TRUE), nrow=W1[1], ncol=W1[2], byrow=TRUE )
    # wgtM <- matrix( 1, nrow=W1[1], ncol=W1[2], byrow=TRUE )
    wgtM <- wgt / wgtM

    wgt <- wgtM

    # missing indicator matrix: 1 - no missings
    missM <- 0.5 * ( (1-is.na(lambda))+ (1- is.na(nu)) )
    wgt <- wgt * missM
    wgt[ missM==0 ] <- 0
    numb_items <- rowSums(missM)

    lambda[ missM==0 ] <- mean( lambda, na.rm=TRUE )
    nu[ missM==0 ] <- mean( nu, na.rm=TRUE )
    group.combis <- t( utils::combn(x=G, m=2))
    group.combis <- rbind( group.combis, group.combis[,c(2,1) ] )
    group_combis <- group.combis

    #--- initial estimates means and SDs
    psi0 <- apply(lambda, 1, stats::median)
    alpha0 <- apply(nu, 1, stats::median, na.rm=TRUE)
    psi0 <- psi0 / psi0[1]
    alpha0 <- alpha0 - alpha0[1]
    if ( ! is.null( psi0.init) ){ psi0 <- psi0.init }
    if ( ! is.null( alpha0.init) ){ alpha0 <- alpha0.init }
    lambda <- as.matrix(lambda)
    wgt <- as.matrix(wgt)

    wgt_combi <- matrix(NA, nrow=nrow(group.combis), ncol=ncol(lambda) )
    for (ii in 1L:I){
        wgt_combi[,ii] <- wgt[ group.combis[,1], ii]*wgt[ group.combis[,2], ii]
    }

    group_combis <- group.combis-1
    G1 <- G-1
    if (overparam){
        G1 <- G
    }
    ind_alpha <- seq_len(G1)
    ind_psi <- G1 + ind_alpha
    if (meth==0){
        ind_alpha <- seq_len(G)
        ind_psi <- G + seq_len(G1)
    }

    #-- define optimization functions
    ia_fct_optim <- function(x, lambda, nu, overparam, eps, meth_,
                        item_wise=FALSE ){
        res <- invariance_alignment_define_parameters(x=x, ind_alpha=ind_alpha,
                        ind_psi=ind_psi, reparam=reparam, meth=meth_)
        alpha0 <- res$alpha0
        psi0 <- res$psi0
        val <- sirt_rcpp_invariance_alignment_opt_fct( nu=nu, lambda=lambda,
                        alpha0=alpha0, psi0=psi0, group_combis=group_combis,
                        wgt=wgt, align_scale=align.scale,
                        align_pow=align.pow, eps=eps, wgt_combi=wgt_combi, type=type,
                        reparam=FALSE, meth=meth_)
        if (item_wise){
            val <- val$fopt_item
        } else {
            val <- val$fopt
        }
        if (overparam | meth==0 ){
            G <- nrow(lambda)
            fac <- sum(wgt[,1]) / 1000
            val <- val+fac*sum(x[1L:G]^2)
        }
        return(val)
    }

    ia_grad_optim <- function(x, lambda, nu, overparam, eps, meth_){
        res <- invariance_alignment_define_parameters(x=x, ind_alpha=ind_alpha,
                            ind_psi=ind_psi, reparam=reparam)
        alpha0 <- res$alpha0
        psi0 <- res$psi0
        grad <- sirt_rcpp_invariance_alignment_opt_grad( nu=nu, lambda=lambda,
                        alpha0=alpha0, psi0=psi0, group_combis=group_combis, wgt=wgt,
                        align_scale=align.scale, align_pow=align.pow, eps=eps,
                        wgt_combi=wgt_combi, type=type, reparam=reparam, meth=meth_)
        grad <- grad[-c(1,G+1)]
        return(grad)
    }
    if ( align.pow[1]==0){
        ia_grad_optim <- NULL
    }
    if ((meth>=2) | num_deriv ){
        ia_grad_optim <- NULL
    }

    #** evaluate optimization function at initial solution
    x0 <- c( alpha0[-1], psi0[-1] )
    if (overparam){
        x0 <- c(alpha0, psi0)
    }
    if (meth==0){
        x0 <- c( alpha0, psi0[-1] )
    }
    fct_optim_inits <- ia_fct_optim(x=x0, lambda=lambda, nu=nu,
                            overparam=overparam, eps=eps, meth=meth )

    #* estimate alignment parameters
    min_val <- .01
    GL <- G1
    par <- c( alpha0[-1], psi0[-1] )
    if (overparam){
        GL <- G
        par <- c( alpha0, psi0 )
    }
    if (meth==0){
        par <- c( alpha0, psi0[-1] )
    }
    lower <- c(rep(-Inf,GL), rep(min_val, GL))
    if (reparam){
        grad_optim <- NULL
    }
    if (meth==0){
        lower <- c(rep(-Inf,G), rep(min_val, GL))
    }
    #* define sequence of epsilon values
    eps_vec <- 10^eps_grid
    eps_vec <- sirt_define_eps_sequence(eps=eps, eps_vec=eps_vec)

    #- optimize (with useful starting values)
    lambda1 <- lambda
    nu1 <- nu
    est_loop <- 1
    psi_list <- list()

    while(est_loop>=1){
        for (eps in eps_vec){
            res_optim <- sirt_optimizer(optimizer=optimizer, par=par, fn=ia_fct_optim,
                                grad=ia_grad_optim, lower=lower, hessian=FALSE,
                                lambda=lambda1, nu=nu1, overparam=overparam,
                                eps=eps, meth_=meth, ...)
            par <- res_optim$par
            res <- invariance_alignment_define_parameters(x=res_optim$par,
                            ind_alpha=ind_alpha, ind_psi=ind_psi, reparam=reparam,
                            meth=meth)
            alpha0 <- res$alpha0
            psi0 <- res$psi0
        } # end eps grid loop
        if (meth%in%c(0,0.5,1,2)){
            est_loop <- 0
        }
        if (meth>=3){
            psi_list[[est_loop]] <- psi0
            if (est_loop==2){
                est_loop <- 0
            }
            if (est_loop==1){
                est_loop <- est_loop + 1
                # nu1 <- nu/sqrt(mod1)*psi_list[[1]]
                #=> likely dead code
            }
        }
    }  # end while loop

    if (meth==3){
        psi0 <- psi_list[[1]]
    }

    #-- standard error computation (if requested)
    TEjkbc <- LEjkbc <- TEjk <- LEjk <- TE <- TEbc <- LEbc <- LE <- SE <- NULL
    vcov0 <- vcov
    V_TEjkbc <- V_TEjk <- V_LEjk <- V_TEbc <- V_TE <- V_LEbc <- V_LE <- NULL
    V_LEjkbc <- NULL

    if (! is.null(vcov)){

        names(par) <- c( paste0('alpha',2L:G), paste0('psi',2L:G) )
        NP <- length(par)

        #- gradient computation
        ia_grad_optim_num <- function(x, lambda, nu, overparam, eps, meth=1,
                                h=1e-4, item_wise=FALSE){
            NP <- length(x)
            par <- x
            if (!item_wise){
                grad <- rep(0,NP)
                names(grad) <- names(x)
            } else {
                grad <- matrix(0,nrow=NP, ncol=I)
                rownames(grad) <- names(x)
            }
            args <- list(x=par, lambda=lambda, nu=nu, overparam=overparam, eps=eps,
                                meth_=meth, item_wise=item_wise)
            for (pp in 1L:NP){
                args$x <- mgsem_add_increment(x=par, h=h, i1=pp)
                f1 <- do.call(what=ia_fct_optim, args=args)
                args$x <- mgsem_add_increment(x=par, h=-h, i1=pp)
                f2 <- do.call(what=ia_fct_optim, args=args)
                der <- (f1-f2)/(2*h)
                if (!item_wise){
                    grad[pp] <- der
                } else {
                    grad[pp,] <- der
                }
            }
            return(grad)
        }

        #- compute hessian with respect to par
        h <- 1e-4
        hess_par <- matrix(NA, nrow=NP, ncol=NP)
        rownames(hess_par) <- colnames(hess_par) <- names(par)
        args <- list(x=par, lambda=lambda, nu=nu, overparam=overparam, eps=eps,
                        meth=meth)
        pp <- 1
        for (pp in 1L:NP){
            args$x <- mgsem_add_increment(x=par, h=h, i1=pp)
            f1 <- do.call( what=ia_grad_optim_num, args=args)
            args$x <- mgsem_add_increment(x=par, h=-h, i1=pp)
            f2 <- do.call( what=ia_grad_optim_num, args=args)
            hess_par[,pp] <- (f1-f2)/(2*h)
        }
        #- compute hessian with respect to lambda and nu
        hess_theta <- matrix(NA, nrow=NP, ncol=2*G*I)
        rownames(hess_theta) <- names(par)
        colnames(hess_theta) <- rownames(vcov)

        pp <- 1

        args <- list(x=par, lambda=lambda, nu=nu, overparam=overparam, eps=eps)

        for (gg in 1L:G){
            for (subs in c('lambda','nu')){
                for (ii in 1L:I){
                    if (subs=='lambda'){
                        z <- lambda
                    } else {
                        z <- nu
                    }
                    args[[subs]] <- mgsem_add_increment(x=z, h=h, i1=gg, i2=ii)
                    f1 <- do.call( what=ia_grad_optim_num, args=args)
                    args[[subs]] <- mgsem_add_increment(x=z, h=-h, i1=gg, i2=ii)
                    f2 <- do.call( what=ia_grad_optim_num, args=args)
                    hess_theta[,pp] <- (f1-f2)/(2*h)
                    pp <- pp+1
                }
            }
        }

        H1 <- MASS::ginv(hess_par)
        A <- H1 %*% hess_theta
        vcov <- A %*% vcov0 %*% t(A)
        rownames(vcov) <- colnames(vcov) <- names(par)
        SE <- sqrt_diag(x=vcov)

        #--- linking error based on M-estimation
        if (le){
            args1 <- list(x=par, lambda=lambda, nu=nu, overparam=overparam, eps=eps,
                            meth=meth, item_wise=TRUE)

            # item-wise gradient
            grad_item <- do.call( what=ia_grad_optim_num, args=args1)
            # variance matrix (meat matrix) for linking error estimation
            M <- 0*vcov
            Mbc <- M
            for (ii in 1L:I){
                gii <- grad_item[,ii]
                Mii <- outer(gii,gii)
                ind <- seq(ii,2*I*G,I)
                vcov0_ii <- vcov0[ind, ind]
                hess_theta_ii <- hess_theta[,ind]
                Mbcii <- hess_theta_ii %*% vcov0_ii %*% t(hess_theta_ii)
                M <- M + Mii
                Mbc <- Mbc + Mbcii
            }
            V_LE <- I/(I-1) * H1 %*% M %*% t(H1)
            V_LEbc <- I/(I-1) * H1 %*% (M-Mbc) %*% t(H1)
            LE <- sqrt_diag(x=V_LE, names=names(par))
            LEbc <- sqrt_diag(x=V_LEbc, names=names(par))

            V_TE <- vcov + V_LE
            V_TEbc <- vcov + V_LEbc
            TE <- sqrt(SE^2 + LE^2)
            TEbc <- sqrt(SE^2 + LEbc^2)

            #--- jackknife linking error

            hess_par_item <- list()
            for (ii in 1L:I){
                hess_par_item[[ii]] <- 0*hess_par
            }

            H0 <- 0*hess_par
            pp <- 1
            for (pp in 1L:NP){
                args1$x <- mgsem_add_increment(x=par, h=h, i1=pp)
                f1 <- do.call( what=ia_grad_optim_num, args=args1)
                args1$x <- mgsem_add_increment(x=par, h=-h, i1=pp)
                f2 <- do.call( what=ia_grad_optim_num, args=args1)
                der <- (f1-f2)/(2*h)
                for (ii in 1L:I){
                    hess_par_item[[ii]][,pp] <- der[,ii]
                }
            }
            for (ii in 1L:I){
                H0 <- H0 + hess_par_item[[ii]]
            }
            # hess_par=H0

            # jk le
            jkfac <- (I-1)/I
            estdiff_jk <- matrix(NA, nrow=I, ncol=NP)
            V_LEjkbc <- 0*V_LE

            for (ii in 1L:I){
                cii <- grad_item[,ii]
                Aii <- hess_par - hess_par_item[[ii]]
                A1ii <- MASS::ginv(Aii)
                parii <- A1ii %*% cii
                estdiff_jk[ii,] <- parii[,1]
                # correction: see LESL paper
                U <- A
                ind <- seq(ii, 2*I*G, I )
                HTii <- hess_theta
                HTii[,ind] <- 0
                Ui <- A1ii %*% HTii
                Wi <- Ui-U
                V_LEjkbc <- V_LEjkbc + jkfac*( Wi %*% vcov0 %*% t(Wi) )
            }
            V_LEjk <- jkfac*crossprod(estdiff_jk)
            LEjk <- sqrt_diag(x=V_LEjk, names=names(par))
            V_LEjkbc <- V_LEjk - V_LEjkbc
            LEjkbc <- sqrt_diag(x=V_LEjkbc, names=names(par))

            TEjk <- sqrt( SE^2 + LEjk^2)
            TEjkbc <- sqrt( SE^2 + LEjkbc^2)

        }  # end le

    } # end vcov

    # center parameters
    res <- invariance_alignment_center_parameters(alpha0=alpha0, psi0=psi0,
                    center=center, meth=meth )
    alpha0 <- res$alpha0
    psi0 <- res$psi0

    # define aligned parameters
    res <- sirt_rcpp_invariance_alignment_opt_fct( nu=nu, lambda=lambda, alpha0=alpha0,
                        psi0=psi0, group_combis=group_combis, wgt=wgt,
                        align_scale=align.scale, align_pow=align.pow, eps=eps,
                        wgt_combi=wgt_combi, type=type, reparam=reparam, meth=1)
    lambda.aligned <- invariance_alignment_process_parameters(par.aligned=res$lambda,
                            par=lambda0)
    nu.aligned <- invariance_alignment_process_parameters(par.aligned=res$nu, par=nu0)
    fopt <- res$fopt

    #**** calculate item statistics and R-squared measures
    # groupwise aligned loading
    # average aligned parameter
    c1 <- invariance_alignment_aligned_parameters_summary(x=lambda.aligned,
                    label='lambda')
    c2 <- invariance_alignment_aligned_parameters_summary(x=nu.aligned, label='nu')
    itempars.aligned <- data.frame(c1, c2, row.names=colnames(lambda) )
    M.lambda_matr <- sirt_matrix2( itempars.aligned$M.lambda, nrow=G)
    M.nu_matr <- sirt_matrix2( itempars.aligned$M.nu, nrow=G)
    lambda.resid <- lambda.aligned - M.lambda_matr
    nu.resid <- nu.aligned - M.nu_matr

    # R-squared measures
    Rsquared.invariance <- c(NA,NA)
    names(Rsquared.invariance) <- c('loadings', 'intercepts' )
    expl <- psi0 * M.lambda_matr
    Rsquared.invariance['loadings'] <- sirt_rsquared(x=lambda, expl=expl)
    expl <- M.nu_matr + alpha0 * M.lambda_matr
    Rsquared.invariance['intercepts'] <- sirt_rsquared(x=nu, expl=expl)

    # correlations aligned parameters
    rbar <- c( invariance_alignment_calc_corr(t(lambda.aligned)),
                    invariance_alignment_calc_corr(t(nu.aligned)) )
    es.invariance <- rbind( Rsquared.invariance, sqrt(1-Rsquared.invariance), rbar)
    rownames(es.invariance) <- c('R2', 'sqrtU2', 'rbar')

    pars <- data.frame(alpha0=alpha0, psi0=psi0)
    rownames(pars) <- rownames(lambda)

    s2 <- Sys.time()
    time_diff <- s2-s1

    #----- OUTPUT:
    res <- list( pars=pars, itempars.aligned=itempars.aligned,
            es.invariance=es.invariance, center=center, lambda.aligned=lambda.aligned,
            lambda.resid=lambda.resid, nu.aligned=nu.aligned, lambda=lambda0,
            nu=nu0, nu.resid=nu.resid, fopt=fopt, align.scale=align.scale,
            align.pow=align.pow0, res_optim=res_optim, eps=eps, wgt=wgt,
            miss_items=missM, numb_items=numb_items, vcov=vcov, V_LE=V_LE,
            V_LEbc=V_LEbc, V_TE=V_TE, V_LEjk=V_LEjk, V_LEjkbc=V_LEjkbc,
            SE=SE, LE=LE, LEbc=LEbc, LEjk=LEjk, LEjkbc=LEjkbc,
            TE=TE, TEbc=TEbc, TEjk=TEjk, TEjkbc=TEjkbc,
            fct_optim_inits=fct_optim_inits, fixed=fixed, meth=meth0,
            meth_internal=meth, s1=s1, s2=s2, time_diff=time_diff, CALL=CALL)
    class(res) <- 'invariance.alignment'
    return(res)
}

