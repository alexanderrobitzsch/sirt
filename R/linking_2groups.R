## File Name: linking_2groups.R
## File Version: 0.139

linking_2groups <- function(pars, method, type="asymm", pow=2, eps=0.001,
            simultan=FALSE, Theta=seq(-6,6,len=101), wgt=NULL, par_init=NULL,
            optimizer="nlminb", control_optimizer=list() )
{
    if (is.null(wgt)){ wgt <- rep(1, length(Theta) ) }

    #- process item parameters
    pars <- pars[, c('a1', 'b1', 'a2', 'b2' ) ]
    pars <- as.matrix(na.omit(pars))
    I <- nrow(pars)

    #- initial parameters
    if (is.null(par_init)){
        par <- c(mu=0,sig=1)
        if (simultan){
            par1 <- linking_2groups_vector_with_names(x=pars[,1],
                            names=paste0('a_',1L:I))
            par2 <- linking_2groups_vector_with_names(x=pars[,2],
                            names=paste0('b_',1L:I))
            par <- c(par, par1, par2)
        }
    } else {
        par <- par_init
    }

    #-- define linking functions
    linking_functions <- list(
        SL=linking_2groups_stocking_lord_fun,
        Hae=linking_2groups_haebara_fun
    )

    grad_functions <- list(
        SL=linking_2groups_stocking_lord_grad,
        Hae=linking_2groups_haebara_grad
    )

    link_fun <- linking_functions[[method]]
    grad_fun <- grad_functions[[method]]  # Will be NULL if method !='Hae'

    #- check gradient
    check_gradient <- FALSE
    # check_gradient <- TRUE

    if (check_gradient){
        requireNamespace('miceadds')
        par[1] <- 1
        par[2] <- 1.2
        if (simultan){
            set.seed(78)
            ind <- 2+I+1L:I
            par2[ind] <- par2[ind] + stats::rnorm(I, sd=1)
        }

        args2 <- list(x=par, pars=pars, Theta=Theta, wgt=wgt, type=type,
                            pow=pow, eps=eps, simultan=simultan )
        f1 <- do.call( what=link_fun, args=args2)
        miceadds::Revalpr('f1')
        h <- 1e-4
        grad <- linking_2groups_numerical_gradient(fun=link_fun, args=args2, h=h)
        miceadds::Revalpr('grad')
        f0 <- do.call( what=grad_fun, args=args2)
        miceadds::Revalpr('f0')
        miceadds::Revalpr('f0-grad')
        miceadds::Revalpr_maxabs('f0','grad')
        stop()
    }

    #* arguments for optimizer
    args <- list( optimizer=optimizer, par=par, fn=link_fun, grad=grad_fun,
                        method='L-BFGS-B', hessian=FALSE, control=control_optimizer,
                        pars=pars, Theta=Theta, wgt=wgt, type=type,
                        pow=pow, eps=eps, simultan=simultan)

    #* optimize link function
    res <- do.call(what=sirt_optimizer, args=args)

    #-- output
    return(res)

}
