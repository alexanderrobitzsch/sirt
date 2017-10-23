## File Name: Q3.R
## File Version: 1.14
## File Last Change: 2017-10-23 15:31:08
 
###############################################################
# Yen's Q3 statistic (1984)
Q3 <- function( dat , theta , b , progress = TRUE )
{
    I <- ncol(dat)
	if (progress){
		cat("Yen's Q3 Statistic based on an estimated theta score \n*** " )		
		cat(I , "Items | " )
	}
    # expected probability
    expected <- .prob.rasch( theta , b )
    # residual 
    residual <- dat - expected  
    # initialize matrix of Q3 values
    I <- ncol(dat)
    q3.matr <- matrix( NA , nrow=I ,ncol= I )
    q3.long <- matrix( NA , nrow=(I-1)*I/2 , ncol=4 )
    colnames(q3.matr) <- rownames(q3.matr) <- colnames(dat)
	q3.matr <- stats::cor( residual , use = "pairwise.complete.obs")		
	nares <- 1 - is.na(residual)
	NIP <- crossprod(nares)
	itempairs <- t( utils::combn( I , 2 ) )
	q3.long[,3] <- q3.matr[ itempairs ]
	q3.long[,4] <- NIP[ itempairs ]
	q3.long[,1] <- colnames(q3.matr)[ itempairs[,1] ]
	q3.long[,2] <- colnames(q3.matr)[ itempairs[,2] ]			
    q3.long <- as.data.frame( q3.long )
    q3.long[,3] <- as.numeric( paste( q3.long[,3] ))
    q3.long <- q3.long[ order( q3.long[,3] ) , ]
    colnames(q3.long) <- c("Item1" , "Item2" , "Q3" , "N" )
    q3.long <-   q3.long[ !is.na( q3.long[,3] ) , ]	
	MQ3 <- mean( q3.long[,3] )
	SDQ3 <- stats::sd( q3.long[,3] )
	Q3.stat <- stats::quantile( q3.long[,3] , prob = c(  .10 , .25 , .50 , .75 , .90  ) )
	Q3.stat <- c("M" = MQ3 , "SD" = SDQ3 , "Min" = min(q3.long[,3] ) , 
						Q3.stat , "Max" = max(q3.long[,3] ) )
	if (progress){
		cat( paste( nrow(q3.long) , "item pairs\n" ) )	
		cat("*** Q3 Descriptives\n")
		print( round(Q3.stat,3))
	}
	#--- OUTPUT
    res <- list( "q3.matrix" = q3.matr , "q3.long" = q3.long ,
				"expected" = expected , "residual" = residual , "Q3.stat" = Q3.stat )
    return(res)
}
#*****************************************************************************************************************

# Q3 <- yen.q3
