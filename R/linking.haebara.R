## File Name: linking.haebara.R
## File Version: 0.423

linking.haebara <- function(itempars, dist="L2", theta=seq(-4,4, length=61),
        optimizer="optim", center=FALSE, eps=1e-5, use_rcpp=TRUE, ...)
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
    is_1pl <- res$is_1pl

    a.orig <- aM
    b.orig <- bM
    prob_theta <- sirt_dnorm_discrete(x=theta, mean=0, sd=1)

    #* define parameter vector
    par <- c(rep(1,NI), rep(0,NI), rep(0,NS-1), rep(1, NS-1))
    parnames <- c( paste0("a_",items), paste0("b_",items), paste0("mu_",studies[-1]),
                    paste0("sigma_",studies[-1]) )
    names(par) <- parnames
    index_a <- 1:NI
    index_b <- NI + 1:NI
    index_mu <- 2*NI + 1:(NS-1)
    index_sigma <- 2*NI + NS - 1 + 1:(NS-1)
    NP <- length(par)

    #-- initial values
    b_mean <- colMeans(bM, na.rm=TRUE)
    par[index_mu] <- - ( b_mean - b_mean[1] )[-1]
    bM_centered <- apply( bM, 2, FUN=function(vv){ vv - mean(vv, na.rm=TRUE) } )
    par[index_b] <- rowMeans(bM_centered, na.rm=TRUE)


    #-- define optimization function
    if (use_rcpp){
        fct_optim_call <- sirt_rcpp_linking_haebara_fct_optim
        grad_optim_call <- sirt_rcpp_linking_haebara_grad_optim
        index_a_ <- index_a - 1
        index_b_ <- index_b - 1
        index_mu_ <- index_mu - 1
        index_sigma_ <- index_sigma - 1
    } else {
        fct_optim_call <- linking_haebara_optim_function_R
        grad_optim_call <- linking_haebara_gradient_function_R
        index_a_ <- index_a
        index_b_ <- index_b
        index_mu_ <- index_mu
        index_sigma_ <- index_sigma
    }

    fct_optim <- function(x){
        a <- x[index_a]
        b <- x[index_b]
        mu <- c(0, x[index_mu])
        sigma <- c(1, x[index_sigma])
        if (is_1pl){
            a <- rep(1,NI)
            sigma <- 1+0*sigma
        }
        args <- list( NI=NI, NS=NS, dist=dist, aM=aM, bM=bM,
                    theta=theta, prob_theta=prob_theta, est_pars=est_pars, wgtM=wgtM,
                    a=a, b=b, mu=mu, sigma=sigma, eps=eps )
        val <- do.call(what=fct_optim_call, args=args)
        return(val)
    }
    grad_optim <- function(x){
        a <- x[index_a]
        b <- x[index_b]
        mu <- c(0, x[index_mu])
        sigma <- c(1, x[index_sigma])
        if (is_1pl){
            a <- rep(1,NI)
            sigma <- 1+0*sigma
        }
        args <- list( NI=NI, NS=NS, dist=dist, aM=aM, bM=bM, theta=theta,
                    prob_theta=prob_theta, est_pars=est_pars, wgtM=wgtM, a=a, b=b,
                    mu=mu, sigma=sigma, eps=eps, index_a=index_a_, index_b=index_b_,
                    index_mu=index_mu_, index_sigma=index_sigma_,
                    parnames=parnames, NP=NP )
        grad <- do.call(what=grad_optim_call, args=args)
        if (is_1pl){
            grad[index_a] <- 0
            grad[index_sigma] <- 0
        }
        names(grad) <- parnames
        return(grad)
    }

    #-- do optimization
    lower <- rep(-Inf, NP)
    lower[index_sigma] <- .001
    res_optim <- sirt_optimizer(optimizer=optimizer, par=par, fn=fct_optim, grad=grad_optim,
                        lower=lower, hessian=FALSE, ... )
    x <- res_optim$par
    a <- x[index_a]
    b <- x[index_b]
    mu <- c(0, x[index_mu])
    sigma <- c(1, x[index_sigma])

    #-- estimate centered parameters
    res <- invariance_alignment_center_parameters(alpha0=mu, psi0=sigma, center=center)
    mu <- res$alpha0
    sigma <- res$psi0
    pars <- data.frame(study=studies, mu=mu, sigma=sigma)
    rownames(pars) <- NULL

    #-- estimated item parameters
    item <- data.frame(item=items, a=a, b=b)
    numb_items <- colSums(est_pars)

    #-- estimated DIF effects
    muM <- sirt_matrix2(x=mu, nrow=NI)
    sigmaM <- sirt_matrix2(x=sigma, nrow=NI)
    a.resid <- aM - a*sigmaM
    # th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)
    b.resid <- bM - ( b - muM )/sigmaM

    #-- output
    s2 <- Sys.time()
    time_diff <- s2-s1
    res <- list( pars=pars, item=item, mu=mu, sigma=sigma, a=a, b=b,
                    a.resid=a.resid, b.resid=b.resid, res_optim=res_optim,
                    est_pars=est_pars, numb_items=numb_items, center=center,
                    dist=dist, eps=eps, NP=NP, use_rcpp=use_rcpp,
                    s1=s1, s2=s2, time_diff=time_diff, CALL=CALL)
    class(res) <- 'linking.haebara'
    return(res)
}
