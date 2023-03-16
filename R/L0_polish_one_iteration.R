## File Name: L0_polish_one_iteration.R
## File Version: 0.354
## File Last Change: 2019-01-14


L0_polish_one_iteration <- function(x, tol, eps=0.01, maxiter=30, type=1)
{
    nr <- nrow(x)
    nc <- ncol(x)
    mod <- NULL
    num0 <- Inf
    elim0 <- NULL
    val0 <- Inf
    ind_x <- ! is.na(x)
    x_update <- x
    for (i in 1L:nr){
        for (j in 1L:nc){
            elim <- c(i,j)
            update <- FALSE
            x1 <- x
            if (ind_x[i,j]){
                x1[elim[1],elim[2]] <- NA
                mod1a <- L1_polish(x=x1, eps=eps, maxiter=maxiter, type=type)
                resid1 <- mod1a$residuals
                num <- sum( abs(resid1) > tol, na.rm=TRUE)
                val_tol <- mean( abs(resid1)*(abs(resid1) <=tol), na.rm=TRUE )
                val <- mean( abs(resid1), na.rm=TRUE )
                if (num < num0){ update <- TRUE }
                if ((num <=num0)&(val<val0)){ update <- TRUE  }
                if (update){
                    num0 <- num
                    val0 <- val
                    elim0 <- elim
                    mod <- mod1a
                    x_update <- x1
                }
            }
        }
    }
    #-- processing
    iterate_further <- TRUE
    max_resid <- max( abs(mod$residuals), na.rm=TRUE)
    if( max_resid <=tol ){
        iterate_further <- FALSE
    }

    #* final least squares estimate
    res1 <- L2_polish(x=x_update)
    row_fin_ls <- res1$row
    col_fin_ls <- res1$col
    wgt <- res1$wgt
    N_elim <- sum(1-wgt)

    #--- output
    res <- list(num=num0, elim=elim0, val=val0, overall=mod$overall,
                    row=mod$row, row_fin_ls=row_fin_ls, col=mod$col,
                    col_fin_ls=col_fin_ls, residuals=mod$residuals,
                    x_update=x_update, tol=tol, iterate_further=iterate_further,
                    max_resid=max_resid, wgt=wgt, N_elim=N_elim)
    return(res)
}


