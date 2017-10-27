## File Name: rm_hrm_est_c_rater.R
## File Version: 0.11

###################################################			
# c.rater
rm_hrm_est_c_rater <- function(  c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=1,theta.k ,
					msteps , mstepconv , est.c.rater , prob.item ,
					c.rater.fixed, c.rater0 )
{
    h <- numdiff.parm
	if (est.c.rater=="r"){ diffindex <- rater.index }
	if (est.c.rater=="i"){ diffindex <- item.index }
	if (est.c.rater=="e"){ diffindex <- rep(1,I) }	
	if (est.c.rater=="a"){ diffindex <- 1:I }		
	RR <- I/VV
	cat("  M steps c.rater parameter    |")
	it <- 0
	conv1 <- 1000
	se.c.rater <- 0 * c.rater
	Q0 <- 0 * c.rater

	#--- input calcprobs
	args <- list( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
						d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, 
						prob.item=prob.item, prob.rater=NULL )	
	
	#--- begin M-steps
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){
		c.rater11 <- b0 <- c.rater
		
		for (kk in 1:K){	
			Q1 <- Q0
			Q1[,kk] <- 1			
			
			args$c.rater <- c.rater11
			pjk <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total	
			args$c.rater <- c.rater11 + h*Q1
			pjk1 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total	
			args$c.rater <- c.rater11 - h*Q1
			pjk2 <- do.call(what=rm_hrm_calcprobs, args=args)$prob.total	
			
			#-- increments
			res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex, 
						max.increment=max.b.increment, numdiff.parm=numdiff.parm ) 								
			c.rater[,kk] <- c.rater[,kk] + res$increment[diffindex]
			se.c.rater[,kk] <- ( sqrt( abs(-1/res$d2) ) )[diffindex ]
			if (kk>1){ 
				ind <- which( c.rater[,kk] < c.rater[,kk-1] )
				if (length(ind)>0){
					l1 <- c.rater[ind,kk-1]
					c.rater[ind, kk-1] <- c.rater[ind,kk]				
					c.rater[ ind , kk ] <- l1
				}
			}
		}
		if ( ! is.null( c.rater.fixed ) ){
		    c.rater[ c.rater.fixed[,1:2] ] <- c.rater.fixed[,3]
		}
		conv1 <- max( abs( c.rater - b0 ) )
		it <- it+1
		cat("-")  
	}
	#---- end M-steps
	
	#- trim increments
	c.rater <- rm_trim_increments_mstep( parm=c.rater, parm0=c.rater0 , max.increment=max.b.increment )
	
	cat(" " , it , "Step(s) \n")	
    res <- list(c.rater = c.rater , se.c.rater = se.c.rater , ll = sum(res$ll0) )
    return(res)
}				

.rm.hrm.est.c.rater <- rm_hrm_est_c_rater

			
