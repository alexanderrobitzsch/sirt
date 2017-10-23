## File Name: rm_hrm_est_tau_item.R
## File Version: 0.11
## File Last Change: 2017-10-03 17:20:48

		
#####################################################################
rm_hrm_est_tau_item <- function( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1  , theta.k ,
				msteps, mstepconv , tau.item.fixed , prob.rater, tau.item0 )
{
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(c.rater)	
	Q0 <- matrix(0,nrow=VV, ncol=K)
	se.tau.item <- Q0
	cat("  M steps tau.item parameter   |")
	it <- 0
	conv1 <- 1000
	tau.item00 <- tau.item0

	#--- input calcprobs
	args <- list( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
					d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, 
					prob.item=NULL, prob.rater=prob.rater )
	
	#--- begin M-steps
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		tau.item11 <- tau.item0 <- tau.item		
		for (kk in 1:K){
			Q1 <- Q0
			Q1[,kk] <- 1
						
			args$tau.item <- tau.item11			
			res <- do.call(what=rm_hrm_calcprobs, args=args)						
			pjk <- res$prob.total	
			prob.item <- res$prob.item

			args$tau.item <- tau.item11 + h*Q1			
			pjk1 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total
			
			args$tau.item <- tau.item11 - h*Q1			
			pjk2 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total
				
			#-- increment
			res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
						max.increment=max.b.increment, numdiff.parm=numdiff.parm ) 
			increment <- Q1 * matrix( res$increment , nrow=VV , ncol=K)	
			tau.item <- tau.item + increment
			se.tau.item[,kk] <- sqrt(abs(-1/res$d2)	)
		}
		conv1 <- max( abs( tau.item - tau.item0 ) )
		it <- it+1
		cat("-") 
		if ( !is.null(tau.item.fixed) ){
			tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
		}
	}
	#---- end M-steps
	
	#- trim increments
	tau.item <- rm_trim_increments_mstep( parm=tau.item, parm0=tau.item00 , max.increment=max.b.increment )	
	
	cat(" " , it , "Step(s) \n")
	res <- list(tau.item = tau.item , se.tau.item = se.tau.item , 
					ll = sum(res$ll0) , prob.item=prob.item )
    return(res)
}

.rm.hrm.est.tau.item <- rm_hrm_est_tau_item
					
					
