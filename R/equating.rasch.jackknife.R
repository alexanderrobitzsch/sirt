## File Name: equating.rasch.jackknife.R
## File Version: 0.14



#--- computation of linking error using Jackknife
equating.rasch.jackknife <- function( pars.data, display=TRUE, se.linkerror=FALSE,
            alpha1=0, alpha2=0 )
{
    pars.data <- as.data.frame( stats::na.omit( pars.data ) )
    itemunits <- unique( pars.data[,1] )
    N.units <- length( itemunits )
    N.items <- nrow( pars.data )
    pars.data[,4] <- paste("I", 1:N.items,sep="")
    # display
    if (display){
        cat( paste( "Jackknife Equating Procedure (Stocking-Lord)\n",
                        N.items, " Items in ", N.units, " Units\n", sep="") )
    }
    # equating without jackknife
    mod1 <- equating.rasch( pars.data[, c( 4, 2) ], pars.data[, c(4, 3) ] )
    res1 <- data.frame( "unit"=itemunits, "shift"=0, "SD"=0, "linkerror"=0)

    # perform jackknife
    for (nn in 1:N.units){
        pars.data1 <- pars.data[ pars.data[,1] !=itemunits[nn], ]
        mod.nn <- equating.rasch( x=pars.data1[, c(4,2) ], y=pars.data1[, c(4,3) ] )
        res1[ nn, "shift" ] <- mod.nn$B.est$Stocking.Lord
        res1[ nn, "SD" ] <- mod.nn$descriptives$SD

        # Jackknife of the linking error
        if (se.linkerror){
            itemunits.nn <- itemunits[ - nn ]
            l1 <- NULL
            for (ii in itemunits.nn){
                pars.data1.ii <- pars.data1[ paste(pars.data1[,1]) !=ii, ]
                mod.ii <- equating.rasch( x=pars.data1.ii[,c(4,2)], y=pars.data1.ii[,c(4,3)],
                                            alpha1=alpha1, alpha2=alpha2)
                l1 <- c(l1, mod.ii$B.est$Stocking.Lord )
            }
            res1[ nn, "linkerror"] <-  sqrt( ( N.units - 2 ) / ( N.units -1 ) * sum( ( l1 - res1[ nn, "shift" ]  )^2 ) )
        }
        # display progress
        if (display){
            cat( paste( nn, " ", sep="" ) )
            utils::flush.console()
            if ( nn%%10==0){ cat("\n") }
        }
    }
    cat("\n")
    linkerror <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,2] - mod1$B.est$Stocking.Lord )^2 ) )
    se.sd <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,3] - mod1$descriptives$SD )^2 ) )
    if (se.linkerror){
        se.linkerror <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,4] - linkerror )^2 ) )
    } else {
        se.linkerror <- NA
    }
    #--- output
    descriptives <- data.frame( N.items=N.items, N.units=N.units, shift=mod1$B.est$Stocking.Lord,
                        SD=mod1$descriptives$SD, linkerror.jackknife=linkerror,
                        SE.SD.jackknife=se.sd, se.linkerror.jackknife=se.linkerror )
    res <- list( pars.data=pars.data, itemunits=itemunits, descriptives=descriptives )
    return(res)
}

