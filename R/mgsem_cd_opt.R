## File Name: mgsem_cd_opt.R
## File Version: 0.170

mgsem_cd_opt <- function(x, opt_fun_args, tol=1e-4, eps_approx=1e-20,
                maxiter=100, h=1e-4, verbose=TRUE, interval_length=0.025,
                method="lqa")
{
    cd_fun_args <- opt_fun_args
    cd_fun_args$eps_approx <- eps_approx
    cd_fun_args$partable_start <- cd_fun_args$partable
    NP <- max(cd_fun_args$partable$index)

    verbose_print <- round( seq( 1,NP, length=min(NP,10)) )

    if (verbose){
        cat('\n----- Start coordinate descent ---- \n')
    }

    iter <- 1
    x00 <- x
    pp <- 1
    maxdiff <- 1

    while ( ( iter < (maxiter+1) ) & ( maxdiff > tol ) ){

        if (verbose){
            cat('Iter.', iter, ' ')
        }

        x0 <- x

        for (pp in 1:NP){
            partable <- cd_fun_args$partable
            x1 <- x[pp]
            interval <- x1 + interval_length*c(-1,1)
            cd_opt_fun <- function(x1)
            {
                x[pp] <- x1
                ll <- mgsem_opt_fun(x=x, opt_fun_args=cd_fun_args)
                return(ll)
            }
            if (method=='exact'){
                res <- stats::optimize(f=cd_opt_fun, interval=interval, tol=tol)
            }
            if (method=='lqa'){
                cda_fun_args <- cd_fun_args
                cda_fun_args$pen_type <- 'none'

                ll0 <- mgsem_opt_fun(x=x, opt_fun_args=cda_fun_args)
                ll1 <- mgsem_opt_fun(x=mgsem_add_increment(x=x, h=h, i1=pp),
                                        opt_fun_args=cda_fun_args)
                ll2 <- mgsem_opt_fun(x=mgsem_add_increment(x=x, h=-h, i1=pp),
                                        opt_fun_args=cda_fun_args)
                a0 <- ll0
                a1 <- (ll1-ll2)/(2*h)
                a2 <- 0.5*(ll1-2*ll0+ll2)/h^2
                x10 <- x1
                cda_opt_fun <- function(x1)
                {
                    x[pp] <- x1
                    res <- mgsem_cda_opt_evaluate_penalties(x=x, opt_fun_args=cd_fun_args)
                    val0 <- -res$pen_all
                    val <- a0+a1*(x1-x10)+a2*(x1-x10)^2+val0
                    return(val)
                }
                res <- stats::optimize(f=cda_opt_fun, interval=interval, tol=tol)
            }

            v1 <- res$minimum
            x[pp] <- v1
            cd_fun_args$partable <- mgsem_coef2partable(coef=x, partable=partable)
            verbose_pp <- verbose & (pp %in% verbose_print)
            if (verbose_pp){
                cat('.')
                utils::flush.console()
            }
        } # end pp
        maxdiff <- max( abs(x-x0) )
        if (verbose){
            cat( ' | Max. diff.=', round(maxdiff, 7), '\n')
        }

        iter <- iter + 1

    }

    value <- cd_opt_fun(x=x1)

    #--- output
    res <- list(par=x, value=value)
    return(res)
}
