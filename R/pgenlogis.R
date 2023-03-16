## File Name: pgenlogis.R
## File Version: 1.11

#-- generalized logistic distribution
pgenlogis <- function( x, alpha1=0, alpha2=0 )
{
    x <- as.vector(x)
    xp <- x[ x >=0 ]
    xn <- x[ x < 0 ]
    indp <- which( x >=0 )
    indn <- which( x < 0 )
    # positive part of x
    if( ( alpha1 > 0 ) ){ y <-  ( ( exp( alpha1 * xp ) - 1  ) / alpha1  ) }
    if( ( alpha1==0 ) ) { y <-  xp  }
    if( ( alpha1 < 0 ) ) { y <- ( - log( 1 - alpha1 * xp )/alpha1 ) }
    # negative part of x
    if( ( alpha2 > 0 ) ){ y1 <-  - ( ( exp( alpha2 * abs(xn) ) - 1  ) / alpha2  ) }
    if( ( alpha2==0 ) ) { y1 <-  xn  }
    if( ( alpha2 < 0 ) ) { y1 <- (  log( 1 - alpha2 * abs(xn) )/alpha2 ) }
    # rearrange entries
    w <- rep( 0, length(x) )
    if ( length(indp) >  0 ){ w[ indp ] <- y }
    if ( length(indn) >  0 ){ w[ indn ] <- y1 }
    # transform to inverse logistic probability
    y <- stats::plogis(w)
    return(y)
}



