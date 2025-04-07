## File Name: linking_2groups_haebara_fun.R
## File Version: 0.167


#- Haebara optimization function
linking_2groups_haebara_fun <- function(x, pars, Theta, wgt, type="asymm",
        pow=2, eps=0.001, simultan=FALSE, deriv=FALSE )
{
    TP <- length(Theta)
    I <- nrow(pars)
    mat <- matrix(NA, nrow=TP, ncol=I)
    mu <- x[1]
    sigma <- x[2]
    val <- 0
    if (deriv){
        NP <- 2
        val <- rep(0,2)
    }
    wgtM <- sirt_matrix2(wgt, nrow=I)

    a1M <- pars[,1]
    b1M <- pars[,2]
    a2M <- pars[,3]
    b2M <- pars[,4]
    if (simultan){
        a <- x[ 2 + 1L:I ]
        b <- x[ 2 + I + 1L:I ]
        ind_a <- 2+1L:I
        ind_b <- 2+I+1L:I
        if (deriv){
            val <- rep(0, 2+2*I)
        }
    }

    if (type=='asymm' | type=='symm'){
        Theta1 <- sirt_matrix2(sigma*Theta+mu, nrow=I)
        Theta2 <- sirt_matrix2(Theta, nrow=I)
        P1 <- stats::plogis( a1M * (Theta1-b1M) )
        P2 <- stats::plogis( a2M * (Theta2-b2M) )
        if (!simultan){ # simultan=FALSE
            z <- P1-P2
            d <- linking_2groups_power_loss(x=z, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){  # no derivative
                val <- val + sum( wgtM*d )
            } else {  # derivative
                fac0 <- P1*(1-P1)*a1M
                val[1] <- val[1] + sum( wgtM*d*fac0)
                val[2] <- val[2] + sum( wgtM*d*fac0*Theta2)
            }
        } else {  # simultan=TRUE
            P0 <- stats::plogis( a * (Theta2-b) )
            z1 <- P1-P0
            z2 <- P2-P0
            d1 <- linking_2groups_power_loss(x=z1, pow=pow, eps=eps, deriv=deriv)
            d2 <- linking_2groups_power_loss(x=z2, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){ # deriv=FALSE
                val <- val + sum( wgtM*(d1+d2) )
            } else {  # deriv=TRUE
                fac0 <- P1*(1-P1)*a1M
                val[1] <- val[1] + sum(wgtM*d1*fac0)
                val[2] <- val[2] + sum(wgtM*d1*fac0*Theta2)
                # derivatives item discriminations
                facP0 <- P0*(1-P0)
                fac <- facP0*(Theta2-b)
                wgtM_d <- wgtM*(d1+d2)
                val[ind_a] <- val[ind_a] - rowSums( wgtM_d*fac )
                # derivatives item difficulties
                fac <- facP0*a
                val[ind_b] <- val[ind_b] + rowSums( wgtM_d*fac )

            }
        }  # end simultan
    }

    if (type=='symm'){
        Theta1 <- sirt_matrix2(Theta, nrow=I)
        Theta2 <- sirt_matrix2(1/sigma*Theta-mu/sigma, nrow=I)
        P1 <- stats::plogis( a1M * (Theta1-b1M) )
        P2 <- stats::plogis( a2M * (Theta2-b2M) )
        if (!simultan){  # simultan=FALSE
            z <- P1-P2
            d <- linking_2groups_power_loss(x=z, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){ # deriv=FALSE
                val <- val + sum( wgtM*d )
            } else {  # deriv=TRUE
                fac0 <- P2*(1-P2)*a2M
                val[1] <- val[1] + sum( wgtM*d*fac0/sigma ) # mu
                val[2] <- val[2] + sum( wgtM*d*fac0*Theta2/sigma ) # sigma
            }  # end deriv
        } else {      # simultan=TRUE
            P0 <- stats::plogis( a * (Theta2-b) )
            z1 <- P1-P0
            d1 <- linking_2groups_power_loss(x=z1, pow=pow, eps=eps, deriv=deriv)
            z2 <- P2-P0
            d2 <- linking_2groups_power_loss(x=z2, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){ # deriv=FALSE
                val <- val + sum( wgtM*(d1+d2) )
            } else {  # deriv=TRUE
                fac1 <- - P0*(1-P0)*a
                fac2 <- P2*(1-P2)*a2M - P0*(1-P0)*a
                val[1] <- val[1] + sum(-wgtM*(d1*fac1+d2*fac2)/sigma)  # mu
                val[2] <- val[2] + sum(-wgtM*(d1*fac1+d2*fac2)*Theta2/sigma)  # sigma
                #- derivatives item discriminations
                facP0 <- P0*(1-P0)
                fac <- facP0*(Theta2-b)
                wgtM_d <- wgtM*(d1+d2)
                val[ind_a] <- val[ind_a] - rowSums( wgtM_d*fac )
                #- derivatives item difficulties
                fac <- facP0*a
                val[ind_b] <- val[ind_b] + rowSums( wgtM_d*fac )

            }
        } # end simultan
    }
    return(val)
}
