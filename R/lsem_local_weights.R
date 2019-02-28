## File Name: lsem_local_weights.R
## File Version: 0.05

lsem_local_weights <- function(data.mod, moderator.grid, h)
{
    eps <- 1E-8
    N <- length(data.mod)
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
    sd.moderator <- stats::sd( data.mod, na.rm=TRUE)
    bw <- h * sd.moderator * N^(-1/5)
    moderator.density <- stats::density( data.mod, from=min(moderator.grid),
                            to=max(moderator.grid ), n=G )$y
    moderator.density <- data.frame( moderator=moderator.grid,
                            wgt=moderator.density / sum(moderator.density) )

    for (gg in 1:G){
        xgg <- moderator.grid[gg]
        wgt <- stats::dnorm( data.mod, mean=xgg, sd=bw ) /
                    stats::dnorm( xgg, mean=xgg, sd=bw )
        weights[,gg] <- ifelse( wgt < eps, eps, wgt )
    }
    #--- output
    res <- list(weights=weights, N=N, G=G, modgrid_index=modgrid_index,
                sd.moderator=sd.moderator, bw=bw, moderator.density=moderator.density)
    return(res)
}

