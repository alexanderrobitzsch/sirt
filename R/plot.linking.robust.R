## File Name: plot.linking.robust.R
## File Version: 0.04


#--- linking.robust S3 plot method
plot.linking.robust <- function( x,  ... )
{
    graphics::par( mfrow=c(2,1))
    KK <- length(x$k.robust)
    graphics::plot( x$k.robust, x$meanpars[1:KK], type="l", xlab="k",
            ylab="Linking constant", main="Linking constant")
    graphics::points( 0, x$meanpars[1], pch=16, col=3, cex=1.4 )
    graphics::points( x$kopt, x$meanpars.kopt, pch=17, col=2, cex=1.4 )
    #****
    graphics::plot( x$k.robust, x$se[1:KK], type="l",
            main=paste0( "Standard error of linking constant (k_opt=", round(x$kopt, 3 ),")" ),
            xlab="k", ylab="Standard error")
    graphics::points( 0, x$se[1], pch=16, col=3, cex=1.4 )
    graphics::points( x$kopt, x$se.kopt, pch=17, col=2, cex=1.4 )
    graphics::par( mfrow=c(1,1) )
}

