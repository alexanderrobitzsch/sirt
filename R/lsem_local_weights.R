## File Name: lsem_local_weights.R
## File Version: 0.196

lsem_local_weights <- function(data.mod, moderator.grid, h,
        sampling_weights=NULL, bw=NULL, kernel="gaussian")
{
    eps <- 1E-8
    N <- length(data.mod)
    no_sampling_weights <- FALSE
    if (is.null(sampling_weights)){
        sampling_weights <- rep(1,N)
        no_sampling_weights <- TRUE
    }
    sampling_weights <- sampling_weights / sum(sampling_weights) * N
    # select nearest neighbor in moderator group for calculating residuals
    G <- length(moderator.grid)
    modgrid_index <- rep(1,N)
    for (gg in 2:G){
        modgrid_index <- ifelse( abs( data.mod - moderator.grid[ modgrid_index ] ) <
                    abs( data.mod - moderator.grid[ gg ] ),
                    modgrid_index, gg )
    }
    # compute weights for every grid point gg
    weights <- matrix( NA, nrow=N, ncol=G )
    sd.moderator <- TAM::weighted_sd(x=data.mod, w=sampling_weights)
    m.moderator <- TAM::weighted_mean(x=data.mod, w=sampling_weights)
    bw_basis <- sd.moderator * N^(-1/5)
    if (is.null(bw)){
        bw <- h * bw_basis
    } else {
        h <- bw / bw_basis
    }
    weights1 <- sampling_weights / sum(sampling_weights)

    moderator.density <- stats::density( data.mod, weights=weights1,
                            from=min(moderator.grid), to=max(moderator.grid ), n=G )$y
    moderator.density <- data.frame( moderator=moderator.grid,
                            wgt=moderator.density / sum(moderator.density) )
    for (gg in 1:G){
        xgg <- moderator.grid[gg]
        wgt <- lsem_kernel_weights(x=data.mod, x0=xgg, bw=bw, kernel=kernel)
        weights[,gg] <- ifelse( wgt < eps, eps, wgt )
    }

    #--- output
    res <- list(weights=weights, N=N, G=G, modgrid_index=modgrid_index,
                m.moderator=m.moderator, sd.moderator=sd.moderator, bw=bw, h=h,
                moderator.density=moderator.density,
                sampling_weights=sampling_weights, no_sampling_weights=no_sampling_weights)
    return(res)
}

