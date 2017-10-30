## File Name: rm_facets_est_a_item.R
## File Version: 0.11


#####################################################
# estimation of slope parameter for items
rm_facets_est_a_item <- function( b.item , b.rater , Qmatrix , tau.item ,
        VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
        n.ik , numdiff.parm=.001 , max.b.increment=1,theta.k , msteps ,
		mstepconv , a.item.center , a.item.fixed , a_lower = .05 , a_upper = 10 )
{
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(b.rater)
	cat("  M steps a.item parameter   |")
	it <- 0
	conv1 <- 1000	
	#--- args calcprobs
	args <- list( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, 
					a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR ) 			
	#---- begin M-steps
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a.item0 <- a.item
		args$a.item <- a.item				
		pjk <- do.call( what=rm_facets_calcprobs, args=args)
		args$a.item <- a.item + h
		pjk1 <- do.call( what=rm_facets_calcprobs, args=args)
		args$a.item <- a.item - h
		pjk2 <- do.call( what=rm_facets_calcprobs, args=args)
		#-- increments
		res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
					max.increment=max.b.increment, numdiff.parm=numdiff.parm ) 				
		a.item <- a.item + res$increment
		#-- bound parameter estimates
		a.item <- rm_squeeze(x=a.item, lower=a_lower, upper=a_upper )
		#-- fixed item parameters
		if ( ! is.null( a.item.fixed) ){		
			ind <- which( ! is.na( a.item.fixed  ) )
		    a.item[ind] <- a.item.fixed[ind]
            res$d2[ind] <- -1E10			
		}
		#-- center item discriminations
		a.item <- rm_center_vector( vec=a.item, center_type=a.item.center, do_log=TRUE)				
		conv1 <- max( abs( a.item - a.item0 ) )
		it <- it+1
		cat("-") 
	}
	cat(" " , it , "Step(s) \n")
	#--- output
    res <- list(a.item = a.item , se.a.item = sqrt( abs(-1/res$d2 )), ll = sum(res$ll0) )
    return(res)
}			
			
.rm.facets.est.a.item <- rm_facets_est_a_item 			
