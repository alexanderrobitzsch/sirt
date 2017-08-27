## File Name: mcmc_WaldTest.R
## File Version: 0.05
## File Last Change: 2017-01-18 11:02:49

##########################################################
# Wald Test for a set of hypotheses
mcmc_WaldTest <- function( mcmcobj , hypotheses ){
	NH <- length(hypotheses)
	n1 <- ncol(mcmcobj)
	mcmcobj <- mcmc_derivedPars( mcmcobj , hypotheses)
	n2 <- ncol(mcmcobj)
	mcmcobj <- mcmcobj[ , seq(n1+1,n2) ]
	v1 <- mcmc_vcov(mcmcobj)
	s1 <- mcmc_summary(mcmcobj)
	c1 <- s1$MAP
	# compute test statistic
	W <- t(c1) %*% solve(v1) %*% c1
	stat <- c( "chi2" = W , "df" = NH)	
	stat["p"] <- 1 - stats::pchisq( W , df = stat["df"])
	res <- list( "hypotheses_summary" = s1 ,
			  "chisq_stat" = stat
					)
	class(res) <- "mcmc_WaldTest"
	return(res)	
		}
##############################################################		
# summary of Wald Test based on MCMC output
summary.mcmc_WaldTest <- function( object , digits = 3 , ... ){
    cat("Wald Test\n")
	W1 <- sprintf( paste0("%." , digits , "f" ) , object$chisq_stat["chi2"] )

	v1 <- paste0("Chi^2 = " ,  W1 , ", df = " , object$chisq_stat["df"])
	v1 <- paste0( v1 , ", p = " , sprintf( paste0("%." , digits , "f" ) , 
					object$chisq_stat["p"] ) )
	cat(v1)

	cat("\n\nSummary Hypotheses\n")	
	obji <- object$hypotheses_summary
	vars <- c("parameter","MAP","SD", "Q2.5", "Q97.5" , "Rhat","SERatio",
					"effSize" )
	obji <- obji[,vars]
	NO <- ncol(obji)
	for (vv in 2:NO ){
		obji[,vv] <- round( obji[,vv] , digits )
		}
	obji[,NO] <- round( obji[,NO] )	
	print(obji)
			}
##################################################################			
