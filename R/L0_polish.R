## File Name: L0_polish.R
## File Version: 0.131


L0_polish <- function(x, tol, conv=0.01, maxiter=30, type=1, verbose=TRUE)
{
    res <- list(x_update=x, iterate_further=TRUE)
    #-- iterate
    while(res$iterate_further){
        res <- L0_polish_one_iteration(x=res$x_update, tol=tol, type=type, eps=conv)
        if (verbose){
            v1 <- paste0('Interactions detected: ', res$N_elim)
            v2 <- paste0(' | Absolute value residual: ', round(res$max_resid,3) )
            cat(v1, v2, '\n')
            utils::flush.console()
        }
    }
    #--- output
    return(res)
}
