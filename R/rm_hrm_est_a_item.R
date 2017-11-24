## File Name: rm_hrm_est_a_item.R
## File Version: 0.11


#########################################################################
rm_hrm_est_a_item <- function( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1,theta.k ,
				msteps, mstepconv , prob.rater , a.item.fixed, a_lower=.05, a_upper=999,
				a_center_type = 2, a.item0, a.prior)
{
	h <- numdiff.parm
	diffindex <- item.index
	RR <- length(c.rater)
	cat("  M steps a.item parameter     |")
	it <- 0
	conv1 <- 1000	
	a.item00 <- a.item0
	
	#--- input calcprobs
	args <- list( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
					d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, 
					prob.item=NULL, prob.rater=prob.rater ) 

	#--- begin M-steps
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a.item0 <- a.item	
		args$a.item <- a.item
		res <- do.call(what=rm_hrm_calcprobs, args=args)
		pjk <- res$prob.total
		prob.item <- res$prob.item
		args$a.item <- a.item + h
		pjk1 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total
		args$a.item <- a.item - h
		pjk2 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total
		#-- increments
		res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
					max.increment=max.b.increment, numdiff.parm=numdiff.parm, prior=a.prior, value=a.item )
		a.item <- a.item + res$increment
		#-- bound parameter estimates
		a.item <- rm_squeeze(x=a.item, lower=a_lower, upper=a_upper )
		#-- center item discriminations		
		if ( is.null( a.item.fixed )){
			a.item <- rm_center_vector( vec=a.item, center_type=a_center_type, do_log=TRUE)						
		}
		if ( ! is.null( a.item.fixed ) ){
			a.item[ a.item.fixed[,1] ] <- a.item.fixed[,2]
		}
		conv1 <- max( abs( a.item - a.item0 ) )
		it <- it+1
		cat("-") 
	}
	#---- end M-steps
	#- trim increments
	a.item <- rm_trim_increments_mstep( parm=a.item, parm0=a.item00 , max.increment=max.b.increment )	

	cat(" " , it , "Step(s) \n")	
	res <- list(a.item = a.item , se.a.item = sqrt( abs(-1/res$d2 )) , 
					ll = sum(res$ll0) , prob.item = prob.item )
	return(res)
}			
###############################################################				


.rm.hrm.est.a.item <- rm_hrm_est_a_item
