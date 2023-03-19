## File Name: lsem_residualize.R
## File Version: 0.433


#**** residualize data
lsem_residualize <- function( data, moderator, moderator.grid,
        lavmodel, h=1.1, bw=NULL, residualize=TRUE, sampling_weights=NULL,
        eps=1E-10, verbose=TRUE, kernel="gaussian", variables_model )
{
    # lavaanify model
    lavaanstr <- sirt_import_lavaan_lavaanify(model=lavmodel)
    vars <- variables_model
    data.mod <- data[, moderator, drop=TRUE ]

    # compute local weights
    res <- lsem_local_weights(data.mod=data.mod, moderator.grid=moderator.grid,
                        h=h, sampling_weights=sampling_weights, bw=bw, kernel=kernel)
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

    sampling_weights <- sampling_weights / sum(sampling_weights) * N
    weights <- weights * sampling_weights

    # residualize
    dat2 <- data
    V <- length(vars)
    residualized_intercepts <- matrix( 0, nrow=G, ncol=V)
    colnames(residualized_intercepts) <- vars
    rownames(residualized_intercepts) <- round( moderator.grid, 3 )
    if (residualize){
        if (verbose){
            cat('** Residualize Data\n')
            utils::flush.console()
        }
        N <- nrow(data)
        for (vv in 1:V){
            var.vv <- vars[vv]
            ind_vv <- which( ! is.na( data[,var.vv] ) )
            y0 <- rep(NA,N)
            for (gg in 1:G){
                x <- dat2[,moderator]
                data1 <- data
                data1$x <- x
                res_formula <- paste0( var.vv, ' ~ x + I(x^2)' )
                weights_gg <- weights[,gg]
                data1$weights_gg <- weights_gg
                mod <- stats::lm( formula=res_formula, data=data1,
                                        weights=weights_gg )
                dfr_pred <- data.frame( x=moderator.grid[gg] )
                m1 <- stats::predict( mod, dfr_pred )
                residualized_intercepts[gg,vv] <- m1
                y <- stats::resid(mod)
                y0[ ind_vv ] <- y
                dat2[, var.vv ] <- ifelse( modgrid_index==gg, y0, dat2[, var.vv ] )
            }
        }
    }

    #--- OUTPUT
    res <- list( resid_vars=vars, data=dat2, weights_grid=weights, bw=bw, h=h,
            moderator.density=moderator.density, sd.moderator=sd.moderator, G=G, N=N,
            residualized_intercepts=residualized_intercepts,
            sampling_weights=sampling_weights, no_sampling_weights=no_sampling_weights,
            m.moderator=m.moderator, residualize=residualize )
    return(res)
}

lsem.residualize <- lsem_residualize
