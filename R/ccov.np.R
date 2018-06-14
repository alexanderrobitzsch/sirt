## File Name: ccov.np.R
## File Version: 1.05



#-----------------------------------------------------
# nonparametric estimation of conditional covariance
ccov.np <- function( data , score , bwscale = 1.1 , thetagrid = seq( -3,3,len=200) ,
        progress=TRUE , scale_score = TRUE ){
            # number of Items I
            I <- ncol(data)
            # z-standardization of score
            if ( scale_score ){
                score <- scale( score )[,1]
            }
            # matrix of item response functions
            if (progress ){
                cat("Pairwise Estimation of Conditional Covariances\n" )
                cat("...........................................................\n" )
                cat("Nonparametric ICC estimation \n " ) }
            icc.items <- matrix( 0 , length(thetagrid) , I )
            if ( I >= 20 ){
                    display <- seq( 1 , I , floor( I/20  ) )[ 2:20 ]
                        } else { display <- 20 }
            i <- 1
            for ( ii in 1:I ){
                x <- score[ ! is.na( data[,ii] )  ]
                y <- data[ ! is.na( data[,ii] ) , ii ]
                icc.items[,ii] <- stats::ksmooth( x  , y , bandwidth = bwscale * length(x)^(-1/5)  ,
                            x.points = thetagrid , kernel="normal")$y
                if ( i < 20 ){ if ( ii == display[i] & progress ){
                            cat( paste( 5*i  , "% " , sep="" ) ) ; i <- i + 1 ;
                                                        if (i == 11){ cat("\n" ) }
                                                        utils::flush.console()} }
                }
            if ( progress){ cat("\n") }
            # weights thetagrid
            wgt.thetagrid <- stats::dnorm(thetagrid)
            wgt.thetagrid <- wgt.thetagrid
            if (progress ){
                cat("...........................................................\n" )
                cat("Nonparametric Estimation of conditional covariances \n " )
                utils::flush.console()
                    }
            # calculation of conditional covariance
            ccov.table <- data.frame( "item1ID" = rep( 1:I , I ) , "item2ID" = rep( 1:I , each = I ) )
            ccov.table <- ccov.table[ ccov.table$item1ID < ccov.table$item2ID , ]
            ccov.table$N <- apply( ccov.table , 1 , FUN = function(ll){
                        sum( rowSums( is.na( data[ , c( ll[1] , ll[2] ) ] ) ) == 0 ) } )
            ccov.table <- ccov.table[ ccov.table$N > 0 , ]
            ccov.table$item1 <- colnames(data)[ ccov.table$item1ID ]
            ccov.table$item2 <- colnames(data)[ ccov.table$item2ID ]
            ccov.table$itempair <- paste( ccov.table$item1 , ccov.table$item2 , sep="-" )
            # smoothing all item pairs
            # calculate conditional covariances
            FF <- nrow( ccov.table )
            ccor.matrix <- ccov.matrix <- prod.matrix <- matrix( 0 , nrow= length(thetagrid ) , ncol = FF )
            ii <- 1
            for (ff in 1:FF){
                if (FF>20){
                    display <- seq( 1 , FF , floor( FF/20  ) )[ 2:20 ]
                } else {
                    display <- seq(1,FF)
                }
                data.ff <- data[ , c( ccov.table[ff,1] , ccov.table[ff,2] ) ]
                which.ff <- which( rowSums( is.na( data.ff ) ) == 0  )
                data.ff <- data.ff[ which.ff , ]
#                y <- data[ ! is.na( data[,ii] ) , ii ]     #       Bug: 2011-07-20
                prod.matrix[,ff] <- stats::ksmooth( x = score[ which.ff]  ,
                                        y = data.ff[,1]*data.ff[,2] ,
                                        bandwidth = bwscale * length(which.ff)^(-1/5)  ,
                                        x.points = thetagrid , kernel="normal")$y
                ccov.matrix[ , ff ] <- prod.matrix[,ff] -  icc.items[, ccov.table[ff,1] ] *
                                        icc.items[, ccov.table[ff,2] ]
                if ( ii < 20 ){
                    if ( ff == display[ii] & progress ){
                        cat( paste( 5*ii  , "% " , sep="" ) )
                        ii <- ii + 1
                        utils::flush.console()
                        if (ii == 11){
                            cat("\n" )
                        }
                    }
                }
            }
            # remove NAs from ccov.matrix
            ccov.matrix[ is.na( ccov.matrix) ] <- 0
            if ( progress ){ cat("\n") }
            # calculate (weighted) conditional covariance
            ccov.table$ccov <- apply( ccov.matrix , 2 , FUN = function(sp){
                        stats::weighted.mean( sp  , wgt.thetagrid ) } )
            res <- list( "ccov.table" = ccov.table , "ccov.matrix" = ccov.matrix ,
                            "data" = data , "score" = score , "icc.items" = icc.items )
            return( res )
            }
#------------------------------------------------------------------------------------------------------------------------------
