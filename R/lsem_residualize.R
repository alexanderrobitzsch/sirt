## File Name: lsem_residualize.R
## File Version: 0.376


#**** residualize data
lsem_residualize <- function( data, moderator, moderator.grid,
        lavmodel, h=1.1, residualize=TRUE, eps=1E-10, verbose=TRUE )
{
    # lavaanify model
    lavaanstr <- lavaan::lavaanify( lavmodel  )
    vars <- unique( c( lavaanstr$rhs, lavaanstr$lhs ) )
    vars <- intersect( colnames(data), vars )
    data.mod <- data[, moderator ]

    # compute local weights
    res <- lsem_local_weights(data.mod=data.mod, moderator.grid=moderator.grid, h=h)
    weights <- res$weights
    modgrid_index <- res$modgrid_index
    N <- res$N
    G <- res$G
    sd.moderator <- res$sd.moderator
    bw <- res$bw
    moderator.density <- res$moderator.density

    # residualize
    dat2 <- data
    V <- length(vars)
    residualized_interceps <- matrix( 0, nrow=G, ncol=V)
    colnames( residualized_interceps ) <- vars
    rownames( residualized_interceps ) <- round( moderator.grid, 3 )

    if (residualize){
        if (verbose){
            cat("** Residualize Data\n")
            utils::flush.console()
        }
        N <- nrow(data)
        for (vv in 1:V){
            # vv <- 1
            var.vv <- vars[vv]
            ind_vv <- which( ! is.na( data[,var.vv] ) )
            y0 <- rep(NA,N)
            for (gg in 1:G){
                # gg <- 1
                x <- dat2[,moderator]
                res_formula <- data[, var.vv ] ~  x + I( x^2 )
                mod <- stats::lm( res_formula, weights=weights[,gg]  )
                m1 <- stats::predict( mod, data.frame( x=moderator.grid[gg] ) )
                residualized_interceps[gg,vv] <- m1
                y <- stats::resid(mod)
                y0[ ind_vv ] <- y
                dat2[, var.vv ] <- ifelse( modgrid_index==gg, y0, dat2[, var.vv ] )
            }
        }
    }
    #--- OUTPUT
    res <- list( resid_vars=vars, data=dat2,
            weights_grid=weights, bw=bw,
            moderator.density=moderator.density,
            sd.moderator=sd.moderator, G=G, N=N,
            residualized_interceps=residualized_interceps )
    return(res)
}

lsem.residualize <- lsem_residualize
