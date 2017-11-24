## File Name: rm_facets_calcprobs.R
## File Version: 0.04


#############################################################################	
# cpp implementation of calculation of facets probabilities
rm_facets_calcprobs <- function( b.item , b.rater , Qmatrix , tau.item ,
		VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
		theta.k , RR )
{
	probs <- rm_facets_calcprobs_cpp( b_item=b.item, b_rater=b.rater, Qmatrix=Qmatrix, tau_item=tau.item, K=K, I=I, 
					TP=TP, a_item=a.item, a_rater=a.rater, item_index=item.index-1, 
					rater_index=rater.index-1, theta_k=theta.k ) 	 
	probs <- array( probs , dim=c(I , K+1 , TP ) )
	return(probs)
}

.rm.facets.calcprobs2 <- rm_facets_calcprobs
