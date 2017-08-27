## File Name: prior_model_pars_CleanString.R
## File Version: 0.02
## File Last Change: 2017-01-18 11:02:51

		
#################################################
# clean string		
prior_model_pars_CleanString <- function( ps ){		
	ps <- gsub( " " , "" , ps )
	ps <- ps[ ps != "" ]
	NP <- length(ps)
	for (pp in 1:NP){		
		# pp <- 1
		ps_pp <- ps[pp]	
		# locate comment symbol
		h1 <- gregexpr(pattern ='#', text= ps_pp)
		if ( h1 > 0){
			ps[pp] <- substring( ps_pp , 1 , h1[[1]] -1)	
		}
	}
	ps <- ps[ ps != "" ]
	return(ps)
}
##################################################	
