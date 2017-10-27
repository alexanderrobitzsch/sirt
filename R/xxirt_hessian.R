## File Name: xxirt_hessian.R
## File Version: 0.23

#############################################
# computation of hessian matrix
xxirt_hessian <- function( object )
{	
	item_list <- object$item_list
	items <- object$items
	Theta <- object$Theta
	ncat <- object$ncat
	partable <- object$partable
	partable_index <- object$partable_index
	dat <- object$dat
	dat_resp <- object$dat_resp
	dat1 <- object$dat1
	resp_index <- object$resp_index
	G <- object$G
	group <- object$group
	maxK <- object$maxK	
	group_index <- object$group_index
	weights <- object$weights
	customTheta <- object$customTheta
	par_items <- object$par_items
	par_Theta <- object$par_Theta
	NPI <- length(par_items)
	NPT <- length(par_Theta)
		
	#***********************
	fct_irt <- function(x){ 
		#*** include free paramaters in partable
		partable <- xxirt_partable_include_freeParameters( partable , x[ 1:NPI ] )		
		#**** include parameter in customTheta
		customTheta$par[ customTheta$est ] <- x[ (NPI+1):(NPI+NPT) ]						
		#*** item probabilities
		probs_items <- xxirt_compute_itemprobs( item_list=item_list , 
								items=items , Theta=Theta , ncat=ncat ,
								partable=partable , partable_index=partable_index )
		#*** compute individual likelihood				
		p.xi.aj <- xxirt_compute_likelihood( probs_items = probs_items , dat=dat , 
							 resp_index=resp_index )														 
		#*** compute prior distribution		
		prior_Theta <- xxirt_compute_priorDistribution( Theta=Theta , 
								  customTheta=customTheta , G=G )					  	  
		#*** compute posterior distribution and expected counts
		res <- xxirt_compute_posterior( prior_Theta=prior_Theta , p.xi.aj=p.xi.aj , 
							group=group ,G=G , weights=weights , dat1=dat1 , 
							dat_resp=dat_resp , maxK=maxK ,group_index = group_index )		
		n.ik <- res$n.ik
		p.aj.xi <- res$p.aj.xi
		N.ik <- res$N.ik
		N.k <- res$N.k
		post_unnorm <- res$post_unnorm							
		dev <- sum( weights * log( rowSums( post_unnorm ) ) )									
		return(dev)
	}
	#*******************************
		
	#--- compute Hessian matrix
	par1 <- xxirt_partable_extract_freeParameters( partable )
	par2 <- xxirt_parTheta_extract_freeParameters( customTheta )
	par <- c(par1 , par2)
	hess <- CDM::numerical_Hessian( par = par , FUN = fct_irt )
	rownames(hess) <- colnames(hess) <- names(par)		
	return(hess)
}
