## File Name: rasch.mml.raschtype.R
## File Version: 2.641



#*******************************************************
# utility function squeeze
squeeze.mml2 <- function( v1, rgvec ){
    v1 <- ifelse( v1 < rgvec[1], rgvec[1], v1 )
    v1 <- ifelse( v1 > rgvec[2], rgvec[2], v1 )
    return(v1)
        }
#*******************************************************


##################################################
# estimation of group means in the 1dim IRT model
.est.mean <- function( dat1.gg, f.yi.qk.gg, X1, pi.k, pi.k0, gg,
                            mean.trait, sd.trait, theta.k, h)
{
    ll0 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k[,gg] ) ) ) )
    pi.k2 <- pi.k0
    pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg]+h, sd=sd.trait[gg]  )
    ll1 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ))
    pi.k2 <- pi.k0
    pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg] -h, sd=sd.trait[gg]  )
    ll2 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    d2[ abs(d2) < 10^(-15) ]  <- 10^(-15)
    d.change <- - d1 / d2
    d.change <- ifelse( abs(d.change) > .1, .1*sign(d.change), d.change )
    return(d.change)
}
##################################################



##################################################
# estimation of group SD's in the 1dim IRT model
.est.sd <- function( dat1.gg, f.yi.qk.gg, X1, pi.k, pi.k0, gg,
                            mean.trait, sd.trait, theta.k, h){
        ll0 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k[,gg] ) ) ) )
        pi.k2 <- pi.k
        pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg], sd=sd.trait[gg] + h )
        ll1 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
        pi.k2 <- pi.k
        pi.k2[,gg] <- sirt_dnorm_discrete( theta.k, mean=mean.trait[gg], sd=sd.trait[gg] - h )
        pi.k2[,gg] <- pi.k2[,gg] / sum( pi.k2[,gg] )
        ll2 <- sum( dat1.gg * log( rowSums( f.yi.qk.gg * outer( X1, pi.k2[,gg] ) ) ) )
        d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
        d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
        d2[ abs(d2) < 10^(-15) ]  <- 10^(-15)
        d.change <- - d1 / d2
        d.change <- ifelse( abs(d.change) > .05, .05*sign(d.change), d.change )
        return(d.change )
                            }
##################################################
