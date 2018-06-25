## File Name: sirt_rbind_fill.R
## File Version: 0.03

## reimplementation of plyr::rbind.fill

sirt_rbind_fill <- function( x, y )
{
    nx <- nrow(x)
    ny <- nrow(y)
    vars <- c( colnames(x), setdiff( colnames(y), colnames(x) ) )
    z <- matrix( NA, nrow=nx+ny, ncol=length(vars) )
    colnames(z) <- vars
    z <- as.data.frame(z)
    z[ 1:nx, colnames(x) ] <- x
    z[ nx + 1:ny, colnames(y) ] <- y
    return(z)
}
