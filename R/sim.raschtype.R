## File Name: sim.raschtype.R
## File Version: 0.23
## File Last Change: 2018-12-30



#---- simulation of Rasch type models
sim.raschtype <- function( theta, b, alpha1=0, alpha2=0, fixed.a=NULL,
        fixed.c=NULL, fixed.d=NULL )
{
    if ( is.null(fixed.a)){ fixed.a <- 1+0*b }
    if ( is.null(fixed.c)){ fixed.c <- 0*b }
    if ( is.null(fixed.d)){ fixed.d <- 1 + 0*b}
    # latent response (subtraction)
    latresp <- TAM::tam_outer( theta, b, op="-" )
    # include slope simulation
    TP <- length(theta)
    latresp <- sirt_matrix2( fixed.a, nrow=TP ) * latresp
    # transformed response
    cM <- sirt_matrix2( fixed.c, nrow=TP )
    dM <- sirt_matrix2( fixed.d, nrow=TP )
    trlat <- pgenlogis( latresp, alpha1=alpha1, alpha2=alpha2 )
    trlat <- cM + ( dM - cM )*trlat
    expprob <- trlat
    N <- nrow(expprob)
    I <- ncol(expprob)
    # define response matrix
    dat.resp <- 1 * ( expprob > matrix( stats::runif(N*I), ncol=I) )
    pot <- max( 2, floor(log10(I) + 1 ) )
    colnames(dat.resp) <- paste( "I", substring( 10^pot + 1:I,2), sep="")
    dat.resp <- as.data.frame(dat.resp)
    return(dat.resp)
}


