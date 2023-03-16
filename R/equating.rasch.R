## File Name: equating.rasch.R
## File Version: 0.244


#---- Equating (linking) in the Rasch model
equating.rasch <- function( x, y, theta=seq( -4, 4, len=100),
        alpha1=0, alpha2=0 )
{
    # Data preparation
    x[,1] <- gsub( " ", "", paste( x[,1] ) )
    y[,1] <- gsub( " ", "", paste( y[,1] ) )
    b.xy <- data.frame( merge( x, y, by.x=1, by.y=1 ) )
    colnames(b.xy) <- c("item", "Itempar.Gr1", "Itempar.Gr2" )
    b.xy <- stats::na.omit( b.xy )
    # mean-mean method
    B.mm <- mean(b.xy[,3]) - mean(b.xy[,2])
    g1 <- .prob.raschtype.genlogis( theta=theta, b=b.xy[,2], alpha1=0, alpha2=0 )
    opt_interval <- 10*c(-1,1)
    #-- Haebara function
    ha <- function(B){
                fct1 <- .prob.raschtype.genlogis( theta=theta, b=b.xy[,2], alpha1=alpha1,
                                    alpha2=alpha2 )
                fct2 <- .prob.raschtype.genlogis( theta=theta, b=b.xy[,3] - B,
                                    alpha1=alpha1, alpha2=alpha2 )
                sum( (fct1 - fct2)^2 )
            }
    B.ha <- stats::optimize( f=ha, interval=opt_interval )$minimum
    # Stocking and Lord Approach
    sl <- function(B){
                fct1 <- .prob.raschtype.genlogis( theta=theta, b=b.xy[,2],
                                alpha1=alpha1, alpha2=alpha2 )
                fct2 <- .prob.raschtype.genlogis( theta=theta, b=b.xy[,3] - B,
                                alpha1=alpha1, alpha2=alpha2 )
                sum( (rowSums( fct1 - fct2 ) )^2 )
            }
    B.sl <- stats::optimize( f=sl, interval=opt_interval )$minimum
    # collect all parameter estimates
    B.est <- data.frame( B.mm, B.ha, B.sl )
    colnames(B.est) <- c("Mean.Mean", "Haebara", "Stocking.Lord")
    # Transformation of item parameters (according to Stocking-Lord)
    b.xy$TransfItempar.Gr1 <- b.xy[,2] + B.est[1,"Stocking.Lord"]
    x[,2] <- x[,2] + B.est[1,"Stocking.Lord"]
    # transformed parameters
    transf.par <- merge( x=x, y=y, by.x=1, by.y=1, all=TRUE )
    colnames(transf.par) <- c("item", "TransfItempar.Gr1", "Itempar.Gr2"  )
    transf.par <- transf.par[ order( paste(transf.par$item ) ), ]
    # calculate variance and linking error
    des <- data.frame( N.Items=nrow(b.xy),
                            SD=stats::sd( b.xy$TransfItempar.Gr1 - b.xy$Itempar.Gr2 ) )
    des$Var <- des$SD^2
    des$linkerror <- sqrt( des["SD"]^2 / des["N.Items"] )[1,1]
    #--- output
    res <- list( "B.est"=B.est, "descriptives"=des,
                        "anchor"=b.xy[, c(1,2,4,3)], "transf.par"=transf.par )
    return(res)
}



