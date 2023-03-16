## File Name: qmc.nodes.R
## File Version: 0.15
## File Last Change: 2019-02-28

qmc.nodes <- function( snodes, ndim )
{
    TAM::require_namespace_msg("sfsmisc")
    r1 <- sfsmisc::QUnif(n=snodes, min=0, max=1, n.min=1, p=ndim, leap=409)
    theta <- as.matrix( stats::qnorm(r1) )
    return(theta)
}
