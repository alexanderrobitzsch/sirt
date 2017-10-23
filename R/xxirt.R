## File Name: xxirt.R
## File Version: 0.903
## File Last Change: 2017-10-02 22:51:38


#########################################################################
# user specified item response model
xxirt <- function( dat , Theta = NULL , itemtype = NULL , customItems = NULL ,
				partable = NULL , customTheta = NULL , group = NULL , weights = NULL ,
				globconv = 1E-6 , conv = 1E-4 , maxit = 200, mstep_iter = 4 , 
				mstep_reltol = 1E-6 , h = 1E-4 , use_grad = TRUE , 
				verbose = TRUE )
{
	#*** preliminaries	
	CALL <- match.call()
	s1 <- Sys.time()
	#*** some data processing of dat		
	res <- xxirt_data_proc(dat=dat , group = group )
    N <- res$N
	G <- res$G
	group <- res$group
	items <- res$items 
    group0 <- res$group0
	groups_unique <- res$groups_unique
	I <- res$I
	maxK <- res$maxK
	ncat <- res$ncat
	weights <- res$weights
	group_index <- res$group_index
	dat_resp <- res$dat_resp
	resp_index <- res$resp_index
	dat1 <- res$dat1 
		
	#*** default Theta
	if ( is.null(Theta) ){
		Theta <- matrix( seq(-6,6,length=21) , ncol=1 )
	}
	TP <- nrow(Theta)
		
	#*** eps - handle numerical instabilities
	eps <- 1E-8		

	# create partable if not provided
	if ( is.null(partable) ){				
		partable <- xxirt_createParTable( dat=dat , itemtype=itemtype , customItems=customItems )				
	}

	# process partable and itemtype				
	res <- xxirt_proc_ParTable( itemtype=itemtype , partable=partable , items=items )
	itemtype <- res$itemtype
	partable <- res$partable
	partable_index <- res$partable_index
	ncat <- res$ncat
	maxK <- res$maxK 
	mstep_method <- res$mstep_method
	item_index <- res$item_index
	dat <- as.matrix(dat) 	

	# create item list												
	item_list <- xxirt_createItemList( customItems = customItems , itemtype=itemtype ,
						items=items, partable = partable )

	# shortcut for calculating expected counts				
	dat1_resp <- xxirt_prepare_response_data(G=G, group_index=group_index, 
						weights=weights, dat1=dat1, dat_resp=dat_resp, maxK=maxK )

	#*** starting values item parameters
	par0 <- xxirt_partable_extract_freeParameters( partable=partable )
	par1 <- xxirt_ThetaDistribution_extract_freeParameters( customTheta=customTheta )
	
	disp <- "...........................................................\n"			
	iter <- 1
	dev <- 1E100
	converged <- FALSE
		
	#*********************************************#
	#************* EM ALGORITHM ******************#
		
	while ( 
			( iter < ( maxit + 1 ) ) & ( ! converged )
				){
				
		if ( verbose){
	       cat(disp)	
	       cat("Iteration" , iter , "   " , paste( Sys.time() ) , "\n" )	
		}			   		
		dev0 <- dev

		#*** item probabilities
		probs_items <- xxirt_compute_itemprobs( item_list=item_list , 
							items=items , Theta=Theta , ncat=ncat ,
							partable=partable , partable_index=partable_index )
							
		#*** compute individual likelihood				
		p.xi.aj <- xxirt_compute_likelihood( probs_items = probs_items , dat=dat , 
							 resp_index=resp_index , dat_resp = dat_resp )
							 
		#*** compute prior distribution		
		prior_Theta <- xxirt_compute_priorDistribution( Theta=Theta , 
							  customTheta=customTheta , G=G )					  
							  
		#*** compute posterior distribution and expected counts
		res <- xxirt_compute_posterior( prior_Theta=prior_Theta , p.xi.aj=p.xi.aj , 
						group=group ,G=G , weights=weights , dat1=dat1 , 
						dat_resp=dat_resp , maxK=maxK ,group_index = group_index,
						dat1_resp=dat1_resp )		
		n.ik <- res$n.ik
		p.aj.xi <- res$p.aj.xi
		N.ik <- res$N.ik
		N.k <- res$N.k
		post_unnorm <- res$post_unnorm
				
		#*** M-step item parameters	    			
		par00 <- par0
		res <- xxirt_mstep_itemParameters( partable=partable, item_list=item_list, 
			          items=items, Theta=Theta , ncat=ncat , 
					  partable_index=partable_index , N.ik=N.ik , 
					  mstep_iter=mstep_iter , par0=par0 , eps = eps,
					  mstep_reltol = mstep_reltol , mstep_method = mstep_method ,
					  item_index = item_index , h = h , use_grad = use_grad )
		ll1 <- res$ll1
		partable <- res$partable
		par0 <- res$par0
			
		#*** M-step theta distribution
		par10 <- par1
		res <- xxirt_mstep_ThetaParameters( customTheta=customTheta, G=G, eps=eps,
				  mstep_iter=mstep_iter  , N.k=N.k , par1=par1 ,
				  mstep_reltol=mstep_reltol , Theta = Theta )		
		ll2 <- res$ll2
		customTheta <- res$customTheta
		par1 <- res$par1	  			
			
		#*** compute deviance
		dev <- - 2 * sum( weights * log( rowSums( post_unnorm ) ) )									
		globconv_temp <- abs( ( - dev + dev0	) / dev0 )					
		
		conv0 <- 0		
		if ( length(par0) > 0){
			conv0 <- max( abs(par0-par00))
		}			
		conv1 <- 0
		if ( !is.null(par1) ){		
			conv1 <- max( abs(par10-par1))	
		}
		conv_temp <- max( conv0 , conv1) 			
        converged <- ( globconv_temp < globconv ) & ( conv_temp < conv )
						
		# print progress	
        if (verbose){			
			cat( paste( "   Deviance = "  , round( dev , 4 ) , 
				if (iter > 1 ){ " | Deviance change = " } else {""} ,
				if( iter > 1){ round( - dev + dev0 , 6 )} else { ""}	,"\n",sep="") )
			cat( paste( "    Maximum item parameter change = " , 
					paste( round( conv0  ,6) , collapse=" " ) , "\n" , sep=""))
			cat( paste( "    Maximum theta distribution parameter change = " , 
					paste( round( conv1  ,6) , collapse=" " ) , "\n" , sep=""))				
			utils::flush.console()						
		}
		iter <- iter + 1
	}
	################### end EM algorithm ####################		
				
	#**** post processing		
				
	#-- parameters
	res <- xxirt_postproc_parameters( partable=partable , customTheta=customTheta , 
		            items=items , probs_items= probs_items )
	par_items <- res$par_items
	par_Theta <- res$par_Theta
	probs_items <- res$probs_items
	par_items_summary <- res$par_items_summary
	par_items_bounds <- res$par_items_bounds
			
	#-- information criteria
	ic <- xxirt_ic( dev=dev , N= sum(weights), par_items=par_items , 
		           par_Theta=par_Theta , I = I , par_items_bounds = par_items_bounds )
		
	#-- compute EAP
	EAP <- xxirt_EAP(p.aj.xi=p.aj.xi , Theta = Theta )
				
		
	#*********************************************#
	# output	
	s2 <- Sys.time()		
	res <- list( partable = partable , par_items = par_items ,
		   			   par_items_summary = par_items_summary , 
					   par_items_bounds = par_items_bounds , 
		               par_Theta = par_Theta , 
					   Theta = Theta, 
					   probs_items = probs_items , 
					   probs_Theta = prior_Theta ,
					   deviance = dev , loglike = -dev/2, ic = ic , item_list=item_list , 
		               customItems = customItems , customTheta = customTheta ,
					   p.xi.aj=p.xi.aj, p.aj.xi=p.aj.xi ,
					   n.ik = n.ik , EAP = EAP ,
					   dat=dat, dat_resp=dat_resp, weights=weights , 
					   item_index = item_index , 
					   G=G , group=group , group_orig = group0 ,
					   ncat=ncat , mstepItem_method = mstep_method ,
					   partable_index = partable_index , items = items , dat1 = dat1 , 
					   group_index = group_index , maxK = maxK , weights = weights ,
					   resp_index = resp_index , 
					   converged = converged , iter = iter-1 , 
					   CALL = CALL , s1 = s1 , s2 = s2
								)								
	class(res) <- "xxirt"
	return(res)		
}
#########################################################################				
