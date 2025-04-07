## File Name: linking_2groups_stocking_lord_fun.R
## File Version: 0.136



#- Stocking-Lord optimization function
linking_2groups_stocking_lord_fun <- function(x, pars, Theta, wgt, type="asymm",
        pow=2, eps=1e-3, simultan=FALSE, deriv=FALSE )
{
    TP <- length(Theta)
    I <- nrow(pars)
    mat <- matrix(NA, nrow=TP, ncol=I)
    mu <- x[1]
    sigma <- x[2]
    a1M <- sirt_matrix2(pars[,1], nrow=TP)
    b1M <- sirt_matrix2(pars[,2], nrow=TP)
    a2M <- sirt_matrix2(pars[,3], nrow=TP)
    b2M <- sirt_matrix2(pars[,4], nrow=TP)

    val <- 0
    if (deriv){
        val <- rep(0,2)
    }

    if (simultan){
        aM <- sirt_matrix2(x[ 2 + 1L:I ], nrow=TP)
        bM <- sirt_matrix2(x[ 2 + I + 1L:I ], nrow=TP)
        ind_a <- 2+1L:I
        ind_b <- 2+I+1L:I
        if (deriv){
            val <- rep(0,2+2*I)
        }
    }

    if (type=='asymm' | type=='symm'){
        Theta1 <- sigma*Theta+mu
        Theta2 <- Theta
        P1 <- stats::plogis( a1M*(Theta1-b1M))
        P2 <- stats::plogis( a2M*(Theta2-b2M))
        rowSums_P1 <- rowSums(P1)
        rowSums_P2 <- rowSums(P2)
        if (! simultan){ # simultan=FALSE
            z <- rowSums_P1 - rowSums_P2
            d <- linking_2groups_power_loss(x=z, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){  # deriv=FALSE
                val <- val + sum( d*wgt )
            } else {  # deriv=TRUE
                fac0 <- P1*(1-P1)*a1M
                val[1] <- val[1] + sum(wgt*d*rowSums(fac0) )  # mu
                val[2] <- val[2] + sum(wgt*d*rowSums(fac0*Theta) )  # sigma
            }
        } else {  # simultan=TRUE
            P0 <- stats::plogis( aM * (Theta2-bM) )
            rowSums_P0 <- rowSums(P0)
            z1 <- rowSums_P1 - rowSums_P0
            d1 <- linking_2groups_power_loss(x=z1, pow=pow, eps=eps, deriv=deriv)
            z2 <- rowSums_P2 - rowSums_P0
            d2 <- linking_2groups_power_loss(x=z2, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){ # deriv=FALSE
                val <- val + sum( (d1+d2)*wgt )
            } else {  # deriv=TRUE
                fac0 <- P1*(1-P1)*a1M
                val[1] <- val[1] + sum(wgt*d1*rowSums(fac0) )  # mu
                val[2] <- val[2] + sum(wgt*d1*rowSums(fac0*Theta) )  # sigma
                #- derivatives item discriminations
                fac1 <- P0*(1-P0)*(Theta2-bM)
                val[ind_a] <- val[ind_a] - colSums(wgt*(d1+d2)*fac1)
                #- derivatives item difficulties
                fac1 <- P0*(1-P0)*aM
                val[ind_b] <- val[ind_b] + colSums(wgt*(d1+d2)*fac1)
            }
        }
    }
    if (type=='symm'){
        Theta1 <- Theta
        Theta2 <- 1/sigma*Theta-mu/sigma
        P1 <- stats::plogis( a1M*(Theta1-b1M))
        P2 <- stats::plogis( a2M*(Theta2-b2M))
        rowSums_P1 <- rowSums(P1)
        rowSums_P2 <- rowSums(P2)
        if (!simultan){ # simultan=FALSE
            z <- rowSums_P1 - rowSums_P2
            d <- linking_2groups_power_loss(x=z, pow=pow, eps=eps, deriv=deriv)
            if (!deriv){ # deriv=FALSE
                val <- val + sum( d*wgt )
            } else {  # deriv=TRUE
                fac0 <- P2*(1-P2)*a2M
                val[1] <- val[1] + sum(wgt*d*rowSums(fac0)/sigma )  # mu
                val[2] <- val[2] + sum(wgt*d*rowSums(fac0*Theta2)/sigma )  # sigma
            }
        } else { # simultan=TRUE
            P0 <- stats::plogis( aM * (Theta2-bM) )
            rowSums_P0 <- rowSums(P0)
            z1 <- rowSums_P1 - rowSums_P0
            d1 <- linking_2groups_power_loss(x=z1, pow=pow, eps=eps, deriv=deriv)
            z2 <- rowSums_P2 - rowSums_P0
            d2 <- linking_2groups_power_loss(x=z2, pow=pow, eps=eps, deriv=deriv)
            wgt_d <- wgt*(d1+d2)
            if (!deriv){  # deriv=FALSE
                val <- val + sum(wgt_d)
            } else {  # deriv=TRUE
                fac1 <- - P0*(1-P0)*aM
                fac2 <- P2*(1-P2)*a2M + fac1
                #* mu
                h1 <- d1*rowSums(fac1)+d2*rowSums(fac2)
                val[1] <- val[1] - sum(wgt*h1/sigma )
                #* sigma
                der <- sum(wgt*h1*Theta2/sigma )
                val[2] <- val[2] - der
                #* derivatives item discriminations
                fac_P0 <- P0*(1-P0)
                fac1 <- fac_P0*(Theta2-bM)
                der <- -colSums(wgt_d*fac1)
                val[ind_a] <- val[ind_a] + der
                #* derivatives item difficulties
                fac1 <- fac_P0*aM
                der <- colSums(wgt_d*fac1)
                val[ind_b] <- val[ind_b] + der
            }
        } # end simultan
    } # end symm

    return(val)
}
