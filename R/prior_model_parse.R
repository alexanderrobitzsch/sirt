## File Name: prior_model_parse.R
## File Version: 0.13
## File Last Change: 2017-01-18 11:02:51

#############################################################
prior_model_parse <- function( prior_model ){
	ps <- strsplit( prior_model , split="\n" , fixed=TRUE)[[1]]
	# clean string
	ps <- prior_model_pars_CleanString( ps )	
	NP <- length(ps)
	prior <- as.list( 1:NP )
	for (pp in 1:NP){
		# pp <- 1
		ps_pp <- ps[pp]
		ps_pp1 <- strsplit( ps_pp , split="~" , fixed=TRUE)[[1]]
		# extract name
		names(prior)[pp] <- ps_pp1[1]
		prior[[pp]] <- as.list(1:2)
		# extract distribution
		ps_pp2 <- ps_pp1[2]
		ps_pp2a <- strsplit( ps_pp2 , split="(" , fixed=TRUE)[[1]]
		prior[[pp]][[1]] <- ps_pp2a[1]
		ps_pp3 <- ps_pp2a[2]
		ps_pp3 <- strsplit( ps_pp3 , split=")" , fixed=TRUE)[[1]][1]	
		ps_pp3 <- strsplit( ps_pp3 , split="," , fixed=TRUE)[[1]]
		NV <- length(ps_pp3)
		prior_pp2 <- as.list( 1:NV )		
		for (vv in 1:NV){
			# vv <- 1
			ps_vv <- ps_pp3[vv]
			h_vv <- strsplit(ps_vv , split="=" , fixed=TRUE)[[1]]
			len_h_vv <- length(h_vv)
			if ( len_h_vv == 1){
				prior_pp2[[vv]] <- suppressWarnings( as.numeric(h_vv) )
			}
			if ( len_h_vv == 2){
				prior_pp2[[vv]] <- suppressWarnings( as.numeric(h_vv[2]) )
				names(prior_pp2)[vv] <- h_vv[1]
			}			
		}
		prior[[pp]][[2]] <- prior_pp2
	}	
	return(prior)
}
###############################################################		
	
