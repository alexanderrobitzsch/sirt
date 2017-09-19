## File Name: mcmc_3pno_testlet_draw_itempars.R
## File Version: 0.02
## File Last Change: 2017-09-19 20:42:04

	
##########################################
# draw item parameters a and b 
mcmc_3pno_testlet_draw_itempars <- function( theta , Z , I , N , weights ,
	gamma.testlet , testletgroups , param , TT , a.testletM)
{
	# define adjusted Z values
	gamma.testletM <- gamma.testlet[ , testletgroups ]
	if (param==1){ Z <- Z - gamma.testletM }
	if (param==3){ Z <- Z - a.testletM*gamma.testletM }	
	if (param==2){ # Z <- Z 
		theta0 <- theta
		Z0 <- Z
	}
	# for parametrization 2, this function must be rewritten
	# because "the theta" is now item specific
	# loop over testlets tt=1,...,TT
	# maybe for TT+1 some adjustment has to be done
	#''''''''''''''''''''''''''''''''''''''''
	# parametrization param=1
	if ( (param==1) | (param==3) ){			
		#--------------
		# sampling without weights
		Xast <- as.matrix( cbind( theta , 1 ) )
		if ( is.null(weights) ){
			Sigma <- solve( crossprod(Xast) )
			# calculate mean
			mj <- Sigma %*% crossprod( Xast , Z )
			mj <- as.matrix( t(mj))		    	
		}
		#--------------
		# sampling with weights						
		if ( ! is.null( weights ) ){
			# compute elements of Xast
			Xast11 <- sum( theta^2 * weights )
			Xast12 <- - sum( theta * weights )
			Xast22 <- sum( weights )
			# compute inverse of Xast
			Xastdet <- Xast11*Xast22 - Xast12^2 
			Xastinv11 <- Xast22 / Xastdet
			Xastinv22 <- Xast11 / Xastdet	
			Xastinv12 <- - Xast12 / Xastdet
			Sigma <- matrix( c(Xastinv11 , Xastinv12 , Xastinv12 , Xastinv22) , 2 ,2 )
			mj <- Sigma %*% crossprod( Xast * weights , Z )
			mj <- as.matrix( t(mj))	
		}		
		#--------------							
		# draw item parameters
		ipars <- sirt_rmvnorm( I , sigma=Sigma ) + mj
		a <- ipars[,1]
		b <- ipars[,2]
	}
	#''''''''''''''''''''''''''''''''''''''''
	# parametrization param=2
	if (param==2){
		a <- rep(NA,I)
		b <- rep(NA,I)
		TTT <- TT
		if ( sum( testletgroups== TT+1 ) > 0 ){
				TTT <- TT + 1 }
		for (tt in 1:TTT){
			#tt <- 1
			theta <- theta0
			Z <- Z0
			ind.tt <- which( testletgroups== tt)
			Itt <- length(ind.tt)		
			theta <- theta0 + gamma.testlet[ , tt]
			Z <- Z[ , ind.tt , drop=FALSE]
			#--------------
			# sampling without weights
			Xast <- as.matrix( cbind( theta , 1 ) )
			if ( is.null(weights) ){
				Sigma <- solve( crossprod(Xast) )
				# calculate mean
				mj <- Sigma %*% crossprod(Xast , Z )
				mj <- as.matrix( t(mj))		    	
								}				
			#--------------
			# sampling with weights						
			if ( ! is.null( weights ) ){
				# compute elements of Xast
				Xast11 <- sum( theta^2 * weights )
				Xast12 <- - sum( theta * weights )
				Xast22 <- sum( weights )
				# compute inverse of Xast
				Xastdet <- Xast11*Xast22 - Xast12^2 
				Xastinv11 <- Xast22 / Xastdet
				Xastinv22 <- Xast11 / Xastdet	
				Xastinv12 <- - Xast12 / Xastdet
				Sigma <- matrix( c(Xastinv11 , Xastinv12 , Xastinv12 , Xastinv22) , 2 ,2 )
				# compute t(Xast) %*% Z (weighted)
				mj <- Sigma %*% crossprod( Xast * weights , Z )
				mj <- as.matrix( t(mj))	
					}		
			#--------------							
			# draw item parameters
			ipars <- sirt_rmvnorm( Itt , sigma=Sigma ) + mj
			a[ind.tt] <- ipars[,1]
			b[ind.tt] <- ipars[,2]
		}			# end testlet tt
	} # end param=2
	#******************************
    res <- list( "a"=a , "b"=b)
	return(res)
}
############################################################


.draw.itempars.3pno.testlet <- mcmc_3pno_testlet_draw_itempars
