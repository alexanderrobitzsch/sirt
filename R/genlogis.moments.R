## File Name: genlogis.moments.R
## File Version: 0.02


#-- moments of generalized logistic distribution
genlogis.moments <- function( alpha1, alpha2)
{
    x0 <- seq(-30, 30, len=30000 )
    y0 <- pgenlogis( x=x0, alpha1=alpha1, alpha2=alpha2 )
    wgt <- y0[-1] - y0[ - length(y0) ]
    wgt <- wgt / sum(wgt)
    out <- ( x0[ -1 ] + x0[ - length(x0) ] ) / 2
    M <- sum( wgt * out )
    SD <- sqrt( sum( wgt*out^2 ) - M^2 )
    moments <- c(M, SD, SD^2 )
    names(moments) <- c('M','SD','Var')
    return(moments)
}
