
#####################################################
# derived parameters for objects of class mcmc
mcmc_derivedPars <- function( mcmcobj , derivedPars ){		
	NP <- length(derivedPars)
	data <- as.data.frame( mcmcobj )
	for (pp in 1:NP){
		# pp <- 1			
		data_pp <- stats::model.matrix( derivedPars[[pp]] , data )
		if (ncol(data_pp) > 1){
			data_pp <- data_pp[,-1]
					}
		data <- as.data.frame( cbind( data , data_pp) )
		colnames(data)[ ncol(data)] <- names(derivedPars)[pp]
					}
    a1 <- attr( mcmcobj , "mcpar")					
	res <- coda::mcmc( 	data= data ,  start = a1[1] , end = a1[2], 
		                thin = a1[3] )
	return(res)
		}
