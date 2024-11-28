## File Name: locpolycor_est_polycor_itempair.R
## File Version: 0.089


locpolycor_est_polycor_itempair <- function(y, ii, jj, data.mod, moderator.grid,
        weights_grid, thresh_ind_list, x0, model="const",
        lower=-.99, upper=.99, EE=99.99, par_init=0.2,
        eps=1e-10, optimizer="optim")
{
    y1 <- y[,ii]
    y2 <- y[,jj]
    data_mod <- data.mod
    G <- length(moderator.grid)
    par0 <- par_init
    if (model=='lin'){
        par0 <- c(par0, 0)
    }
    K1 <- max(y1)
    K2 <- max(y2)
    N <- length(y1)
    x1 <- data_mod

    W <- polycor <- rep(NA,G)

    for (gg in 1L:G){
        x0 <- moderator.grid[gg]
        w <- weights_grid[,gg]
        thresh_ii <- thresh_ind_list[[ii]]
        thresh_jj <- thresh_ind_list[[jj]]
        thresh_ii <- cbind(-EE, thresh_ii, EE)
        thresh_jj <- cbind(-EE, thresh_jj, EE)
        thresh_ii1 <- thresh_ii[ cbind(1L:N, y1+2) ]
        thresh_ii2 <- thresh_ii[ cbind(1L:N, y1+1) ]
        thresh_jj1 <- thresh_jj[ cbind(1L:N, y2+2) ]
        thresh_jj2 <- thresh_jj[ cbind(1L:N, y2+1) ]

        res_optim <- sirt_optimizer(optimizer=optimizer, par=par0,
                            fn=locpolycor_est_polycor_opt_fun,
                            w=w, thresh_ii1=thresh_ii1,
                            thresh_ii2=thresh_ii2, thresh_jj1=thresh_jj1,
                            thresh_jj2=thresh_jj2, eps=eps,
                            x1=data.mod, x0=x0, model=model, method='L-BFGS-B',
                            lower=lower, upper=upper, hessian=FALSE,
                            package='pbivnorm')
        par0 <- res_optim$par
        polycor[gg] <- par0[1]
        W[gg] <- sum(w)
    }

    #-- output
    res <- list( polycor=polycor, ii=ii, jj=jj, W=W, N=N)
    return(res)
}
