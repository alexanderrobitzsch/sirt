## File Name: lsem_residualize.R
## File Version: 0.478


#**** residualize data
lsem_residualize <- function( data, moderator, moderator.grid,
        lavmodel, h=1.1, bw=NULL, residualize=TRUE, sampling_weights=NULL,
        eps=1E-10, verbose=TRUE, kernel="gaussian", variables_model,
        is_imputed=FALSE, Nimp=0 )
{
    # lavaanify model
    lavaanstr <- sirt_import_lavaan_lavaanify(model=lavmodel)
    vars <- variables_model

    # values of moderator variable
    if (! is_imputed){
        data.mod <- data[, moderator, drop=TRUE ]
        Nimp <- 1
    } else {  # use first dataset for imputed data
        data.mod <- ( data[[1]] )[, moderator, drop=TRUE ]
    }

    # compute local weights
    res <- lsem_local_weights(data.mod=data.mod, moderator.grid=moderator.grid,
                        h=h, sampling_weights=sampling_weights, bw=bw, kernel=kernel,
                        is_imputed=is_imputed, Nimp=Nimp, data=data,
                        moderator=moderator)
    weights <- res$weights
    modgrid_index <- res$modgrid_index
    N <- res$N
    G <- res$G
    m.moderator <- res$m.moderator
    sd.moderator <- res$sd.moderator
    bw <- res$bw
    h <- res$h
    moderator.density <- res$moderator.density
    sampling_weights <- res$sampling_weights
    no_sampling_weights <- res$no_sampling_weights

    res0 <- as.list(1L:Nimp)
    data0 <- data

    # residualize
    for (ii in 1L:Nimp){
        dat2 <- data
        if (is_imputed){
            dat2 <- data <- data0[[ii]]
        }

        V <- length(vars)
        residualized_intercepts <- matrix( 0, nrow=G, ncol=V)
        colnames(residualized_intercepts) <- vars
        rownames(residualized_intercepts) <- round( moderator.grid, 3 )
        if (residualize){
            if (verbose & ii==1){
                cat('** Residualize Data\n')
                utils::flush.console()
            }
            N <- nrow(data)
            for (vv in 1L:V){
                var.vv <- vars[vv]
                ind_vv <- which( ! is.na( data[,var.vv] ) )
                y0 <- rep(NA,N)
                for (gg in 1L:G){
                    x <- dat2[,moderator]
                    data1 <- data
                    data1$x <- x
                    res_formula <- paste0( var.vv, ' ~ x + I(x^2)' )
                    if (!is_imputed){
                        weights_gg <- weights[,gg]
                    } else {
                        weights_gg <- (weights[[ii]])[,gg, drop=TRUE]
                    }
                    data1$weights_gg <- weights_gg
                    mod <- stats::lm( formula=res_formula, data=data1,
                                            weights=weights_gg )
                    dfr_pred <- data.frame( x=moderator.grid[gg] )
                    m1 <- stats::predict( mod, dfr_pred )
                    residualized_intercepts[gg,vv] <- m1
                    y <- stats::resid(mod)
                    y0[ ind_vv ] <- y
                    if (!is_imputed){
                        modgrid_index1 <- modgrid_index
                    } else {
                        modgrid_index1 <- modgrid_index[[ii]]
                    }

                    dat2[, var.vv] <- ifelse( modgrid_index1==gg, y0, dat2[, var.vv] )
                }
            }  # end vv
        }  # end residualize=TRUE
        res <- list( resid_vars=vars, data=dat2, weights_grid=weights, bw=bw, h=h,
                moderator.density=moderator.density, sd.moderator=sd.moderator, G=G, N=N,
                residualized_intercepts=residualized_intercepts,
                sampling_weights=sampling_weights,
                no_sampling_weights=no_sampling_weights,
                m.moderator=m.moderator, residualize=residualize,
                is_imputed=is_imputed, Nimp=Nimp)
        res0[[ii]] <- res

    }  # end imputation loop ii

    #--- process output
    if (! is_imputed){
        res <- res0[[1]]
    } else {
        res <- res0[[1]]
        entries <- c('data','residualized_intercepts')
        for (ee in entries){
            v1 <- list()
            for (ii in 1L:Nimp){
                v1[[ii]] <- res0[[ii]][[ee]]
            }
            res[[ee]] <- v1
        }

        #- aggregate residualized intercepts
        ee <- 'residualized_intercepts'
        res[[ee]] <- lsem_aggregate_statistics(x=res[[ee]])

    }  # end process output

    #-- out
    return(res)
}

lsem.residualize <- lsem_residualize
