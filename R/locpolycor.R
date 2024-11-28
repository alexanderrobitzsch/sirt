## File Name: locpolycor.R
## File Version: 0.251


locpolycor <- function(y, data.mod, moderator.grid, h=1.1,
                    model_thresh, model_polycor, kernel="gaussian",
                    eps=1e-10)
{
    #- compute weights
    llw <- lsem_local_weights(data.mod=x, moderator.grid=moderator.grid, h=h,
                kernel=kernel)
    weights_grid <- llw$weights

    #- estimate thresholds
    y <- as.matrix(y)
    I <- ncol(y)
    G <- length(moderator.grid)

    #- item-wise optimization
    thresh_list <- as.list(1L:I)
    K_list <- as.list(1L:I)
    thresh_ind_list <- as.list(1L:I)
    items <- colnames(y)
    names(thresh_list) <- items
    names(thresh_ind_list) <- items
    names(K_list) <- items
    Kmax <- max(y, na.rm=TRUE)
    thresh_stat <- list()
    thresh0 <- matrix(NA, nrow=Kmax, ncol=I)
    rownames(thresh0) <- paste0('t',1L:Kmax)
    colnames(thresh0) <- items
    for (gg in 1L:G){
        thresh_stat[[gg]] <- thresh0
    }

    for (ii in 1L:I){
        cat(paste0( '-- compute thresholds for item ', ii, '\n') )
        utils::flush.console()
        y1 <- y[,ii]
        ind <- which( ! is.na(y1) )
        K <- max(y1)
        thresh_item <- matrix(NA, nrow=K, ncol=G)
        thresh_ind_cbind <- NULL
        data_mod <- data.mod[ind]
        res_item <- list()
        par_init <- NULL
        for (gg in 1L:G){
            x0 <- moderator.grid[gg]
            w <- weights_grid[ind, gg]
            res <- locpolycor_est_thresh_item(y=y1, data.mod=data_mod, x0=x0,
                        w=w, model=model_thresh, par_init=par_init, eps=eps)
            par_init <- res$res_optim$par
            thresh_ind_cbind <- cbind( thresh_ind_cbind, res$thresh_ind)
            thresh_item[,gg] <- res$thresh
            res_item[[gg]] <- res
            thresh_stat[[gg]][1L:K,ii] <- res$thresh
        }

        thresh_list[[ii]] <- thresh_item

        K_list[[ii]] <- K

        #* individual predictions of thresholds
        thresh_ind <- matrix(NA, nrow=N, ncol=K)
        distmat <- abs( outer( data_mod, moderator.grid, '-' ) )
        smallest <- rowKSmallest2.sirt(matr=distmat, K=2)
        smallest <- smallest$smallind[,1L:2]

        # define predictions
        thresh_ind <- res_item[[1]]$thresh_ind
        ind <- which(data_mod > moderator.grid[G])
        if (length(ind)>0){
            thresh_ind[ind, ] <- ( res_item[[G]]$thresh_ind )[ ind,, drop=FALSE ]
        }

        ind <- which( ( data_mod <=moderator.grid[G]) & (data_mod >=moderator.grid[1]))
        sm1 <- smallest[,1]
        sm2 <- smallest[,2]
        mg1 <- moderator.grid[sm1]
        mg2 <- moderator.grid[sm2]
        t3 <- t2 <- t1 <- matrix(NA, nrow=N, ncol=K)
        for (kk in 1L:K){
            mat1 <- cbind( 1L:N, kk+K*(sm1-1) )
            t1[,kk] <- thresh_ind_cbind[ mat1 ]
            mat2 <- cbind( 1L:N, kk+K*(sm2-1) )
            t2[,kk] <- thresh_ind_cbind[ mat2 ]
            # linear interpolation of individual thresholds
            t3[,kk] <- t1[,kk] + ( t2[,kk]-t1[,kk]) * (data_mod-mg1) / (mg2-mg1)
        }
        thresh_ind[ind,] <- t3[ind,]
        thresh_ind_list[[ii]] <- thresh_ind

    }  # end ii

    #-- local polychoric correlations
    polycor0 <- diag(1,I)
    rownames(polycor0) <- colnames(polycor0) <- items
    polycor_stat <- list()
    for (gg in 1L:G){
        polycor_stat[[gg]] <- polycor0
    }
    for (ii in 1L:(I-1) ){
        for (jj in (ii+1):I){
            cat(paste0( '** compute polychoric correlation for item pair (',
                                ii, ',', jj, ')\n') )
            res <- locpolycor_est_polycor_itempair(y=y, ii=ii, jj=jj, data.mod=data.mod,
                            moderator.grid=moderator.grid, weights_grid=weights_grid,
                            model=model_polycor, thresh_ind_list=thresh_ind_list,
                            x0=x0, eps=eps)
            polycor1 <- res$polycor
            for (gg in 1L:G){
                polycor_stat[[gg]][ii,jj] <- polycor1[[gg]]
                polycor_stat[[gg]][jj,ii] <- polycor1[[gg]]
            }
        }
    }

    #-- output
    res <- list(thresh_list=thresh_list, thresh_stat=thresh_stat,
                    polycor_stat=polycor_stat, thresh_ind_list=thresh_ind_list,
                    K_list=K_list, I=I, moderator.grid=moderator.grid,
                    weights_grid=weights_grid)
    return(res)
}
