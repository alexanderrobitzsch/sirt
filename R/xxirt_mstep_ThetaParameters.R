## File Name: xxirt_mstep_ThetaParameters.R
## File Version: 0.06
## File Last Change: 2017-01-18 11:02:56

#######################################################################
xxirt_mstep_ThetaParameters <- function( customTheta , G , eps ,
			mstep_iter , N.k , par1 , mstep_reltol , Theta ){
		like_Theta <- function( x , ... ){
			par1 <- customTheta$par 
			par1[ customTheta$est ] <- x		
			arg_list <- list( "par" = par1 , "Theta" = Theta , "G" = G )
			mod1 <- do.call( customTheta$P , arg_list )
			ll2 <- - sum( N.k * log( mod1 + eps ) )
			NP <- length(customTheta$prior)
			pen <- 0	
			if ( NP > 0 ){			
				for (pp in 1:NP){
					# pp <- 1
					if ( ! is.na( customTheta$prior[pp] ) ) {									
						prior_pp <- customTheta$prior[pp]
						val <- par1[ names(customTheta$prior) ]	
						arg_list <- list( val ,customTheta$prior_par1[pp] , 
												customTheta$prior_par2[pp]  )
						prior_val <- do.call( prior_pp , arg_list )										
						pen <- pen - log( prior_val + eps )
					}  # end if prior
				}  # end pp			
			}
			ll2 <- ll2 + pen
			return(2*ll2)
		}
		#----- end definition likelihood function
					
		# mstep_method <- "L-BFGS-B"		
		mstep_method <- "BFGS"
		
		arg_control <- list(maxit=mstep_iter)		
		if ( mstep_method == "BFGS"){
			arg_control$reltol <- mstep_reltol		
		}
        # method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'		
		arg_list <- list( par = par1 , fn = like_Theta , method = mstep_method ,
				                control=arg_control)
		mod <- do.call( stats::optim , arg_list )						
		par1 <- mod$par
		ll2 <- mod$value						
		customTheta$par[ customTheta$est ] <- par1				
		par1 <- xxirt_ThetaDistribution_extract_freeParameters( customTheta=customTheta )
		res <- list( par1 = par1 , ll2 = ll2, customTheta = customTheta)
		return(res)
}
#######################################################################		
