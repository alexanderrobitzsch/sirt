## File Name: rm.facets.R
## File Version: 4.57

#################################################################
# Facets Model for Raters:
# MML estimation
rm.facets <- function( dat , pid=NULL , rater=NULL ,
	Qmatrix=NULL , theta.k=seq(-9,9,len=30) , 
	est.b.rater=TRUE , est.a.item=FALSE , est.a.rater=FALSE ,
	est.mean = FALSE , tau.item.fixed=NULL , a.item.fixed=NULL , b.rater.fixed=NULL , a.rater.fixed=NULL , 
	b.rater.center = 2, a.rater.center=2, a.item.center=2, 
	max.b.increment=1 , numdiff.parm=.00001 , maxdevchange=.10 ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001, PEM=FALSE, PEM_itermax=maxiter)
{

	CALL <- match.call()
	s1 <- Sys.time()
	
	skillspace <- "normal"
	dat <- as.matrix(dat)	
	if ( is.null(rater)){	
		rater <- rep(1,nrow(dat)) 
		est.b.rater <- FALSE
		est.a.rater <- FALSE		
	}
	if ( is.null(pid)){  
		pid <- seq(1,nrow(dat) ) 
	}
	pcm.param <- FALSE	
	theta.k0 <- theta.k
	pi.k <- sirt_dnorm_discrete(x=theta.k, mean=0, sd=1)	
	
	# process data
	procdata <- res <- rm_proc( dat=dat , rater=rater , pid=pid )
	dat2 <- as.matrix(res$dat2)
	dat2.resp <- as.matrix(res$dat2.resp)
	rater.index1 <- res$rater.index
	dataproc.vars <- res$dataproc.vars
	VV <- res$VV
	RR <- res$RR
	item.index <- res$dataproc.vars$item.index 
	rater.index <- res$dataproc.vars$rater.index
	dat2.ind.resp <- res$dat2.ind.resp	
	
	deviance.history <- rep(NA, maxiter)
	
	# maximum categories
	maxK <- sirt_colMaxs(x=dat)
	
	K <- max( maxK )
	if ( is.null(Qmatrix) ){
		Qmatrix <- matrix( 1:K , nrow=VV , ncol=K , byrow=TRUE)
	}
	TP <- length(theta.k)
	I <- VV*RR

	# center parameters
	if ( ! is.null( b.rater.fixed) ){
	    b.rater.center <- 0		
	}
	if ( ! is.null( a.rater.fixed) ){
	    a.rater.center <- 0		
	}
	if ( ! is.null( a.item.fixed) ){
	    a.item.center <- 0		
	}				
	
	if ( skillspace == "loglinear" ){
		est.mean <- TRUE
	}
	
	# define constraints on tau.item parameters
	# if not all categories are observed
	tau.item.fixed_val <- tau.item.fixed
	tau.item.fixed <- NULL
	if ( min(maxK) < K ){
		tau.item.fixed <-  rm_determine_fixed_tau_parameters( K=K, maxK=maxK, VV=VV )
	}				
					
	# starting values for item difficulties
	b.item <- - stats::qlogis( colMeans( dat , na.rm=TRUE ) / maxK  )
	if ( ! pcm.param ){ 
		b.item <- 0*b.item	
	}
	
	tau.item <- matrix( 0 , nrow=VV , ncol=K )
	rownames(tau.item) <- colnames(dat)
	tau.item <- matrix( seq( -2 , 2 , len=K ) , nrow=VV , ncol=K , byrow=TRUE )

	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
	N <- stats::aggregate( N , list( rater.index ) , sum )[,2]
	M1 <- stats::aggregate( M1 , list( rater.index ) , mean )[,2]		
	b.rater <- - stats::qlogis( M1 / K )
	b.rater <- b.rater - mean( b.rater )
	a.item <- rep(1,VV)
	a.rater <- rep(1,RR)
	if ( ! est.b.rater ){ 
		b.rater <- rep(0,RR) 
	}

	# init standard errors
	se.b.rater <- NA*b.rater
	se.a.rater <- NA*a.rater
	se.a.item <- NA*a.item
	
	#-- preliminaries PEM acceleration
	if (PEM){
		res <- rm_facets_pem_inits( tau.item=tau.item, a.item=a.item, a.rater=a.rater, 
						b.rater=b.rater, skillspace=skillspace, PEM=PEM, a.item.fixed=a.item.fixed, 
						est.a.item=est.a.item ) 
		PEM <- res$PEM
		pem_pars <- res$pem_pars
	#	center_log_a <- res$center_log_a
		pem_parameter_index <- res$pem_parameter_index
		pem_parameter_sequence <- res$pem_parameter_sequence
	}
	
	# inits
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	sigma <- 1
	disp <- "...........................................................\n"	

	b.rater.incr <- max.b.increment
	tau.item.incr  <- max.b.increment
	
active <- TRUE
active <- FALSE
	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter )	){
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
		
		# previous values
		b.item0 <- b.item
		b.rater0 <- b.rater
		tau.item0 <- tau.item
		dev0 <- dev
		sigma0 <- sigma
		a.item0 <- a.item
		a.rater0 <- a.rater
zz0 <- Sys.time()
		# calculate probabilities
		probs <- rm_facets_calcprobs( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, 
						I=I, TP=TP, a.item=a.item, a.rater=a.rater, item.index=item.index, 
						rater.index=rater.index, theta.k=theta.k, RR=RR )					
zz0 <- sirtcat( "  *** facets.calcprob   " , zz0 , active )
		# calculate posterior
		res <- rm_posterior( dat2=dat2, dat2.resp=dat2.resp, TP=TP, pi.k=pi.k, K=K, I=I, probs=probs, dat2.ind.resp=dat2.ind.resp ) 
zz0 <- sirtcat( "  *** posterior   " , zz0 , active )
		f.yi.qk <- res$f.yi.qk
		f.qk.yi <- res$f.qk.yi
		n.ik <- res$n.ik
		N.ik <- res$N.ik
		pi.k <- res$pi.k
		ll <- res$ll
		
		#--- estimate b.rater parameter
		if( est.b.rater ){
			res <- rm_facets_est_b_rater( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, 
							I=I, TP=TP, a.item=a.item, a.rater=a.rater, item.index=item.index, 
							rater.index=rater.index, n.ik=n.ik, numdiff.parm=numdiff.parm, 
							max.b.increment=b.rater.incr, theta.k=theta.k, msteps=msteps, 
							mstepconv=mstepconv, b.rater.center=b.rater.center, 
							b.rater.fixed=b.rater.fixed ) 
			b.rater <- res$b.rater
			se.b.rater <- res$se.b.rater
			b.rater.incr <- abs( b.rater0 - b.rater )
		}
zz0 <- sirtcat( "  *** est b.rater   " , zz0 , active )
												
		#--- estimate tau.item parameters
		res <- rm_facets_est_tau_item( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, 
					I=I, TP=TP, a.item=a.item, a.rater=a.rater, item.index=item.index, 
					rater.index=rater.index, n.ik=n.ik, numdiff.parm=numdiff.parm, 
					max.b.increment=tau.item.incr, theta.k=theta.k, msteps=msteps, 
					mstepconv=mstepconv, tau.item.fixed=tau.item.fixed, 
					tau.item.fixed_val=tau.item.fixed_val ) 
		tau.item <- res$tau.item
		se.tau.item <- res$se.tau.item
		tau.item.incr  <- abs( tau.item0 - tau.item )
zz0 <- sirtcat( "  *** est tau.item   " , zz0 , active )									
		
		#--- estimate a.item parameter
		if (est.a.item){
			res <- rm_facets_est_a_item( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, 
						a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index, 
						n.ik=n.ik, numdiff.parm=numdiff.parm, max.b.increment=1, theta.k=theta.k, msteps=msteps, 
						mstepconv=mstepconv, a.item.center=a.item.center, a.item.fixed=a.item.fixed ) 	
			a.item <- res$a.item
			se.a.item <- res$se.a.item
		}
zz0 <- sirtcat( "  *** est a.item   " , zz0 , active )									
														
		#--- estimate a.rater parameter
		if (est.a.rater){
			res <- rm_facets_est_a_rater( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, 
						a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index, 
						n.ik=n.ik, numdiff.parm=numdiff.parm, max.b.increment=1, theta.k=theta.k, msteps=msteps, 
						mstepconv=mstepconv, a.rater.center=a.rater.center, a.rater.fixed=a.rater.fixed ) 
			a.rater <- res$a.rater
			se.a.rater <- res$se.a.rater
		}
zz0 <- sirtcat( "  *** est a.rater   " , zz0 , active )

		utils::flush.console()	

		#-- update distribution
		res <- rm_smooth_distribution( theta.k=theta.k, pi.k=pi.k, est.mean=est.mean, skillspace=skillspace ) 
		pi.k <- res$pi.k
		mu <- res$mu
		sigma <- res$sigma

		#-- PEM acceleration
		if (PEM){
			res <- rm_facets_pem_acceleration( iter=iter, pem_parameter_index=pem_parameter_index, 
						pem_parameter_sequence=pem_parameter_sequence, a.rater=a.rater, Qmatrix=Qmatrix, 
						tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, b.item=b.item, b.rater=b.rater, 
						item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, dat2=dat2, 
						dat2.resp=dat2.resp, pi.k=pi.k, dat2.ind.resp=dat2.ind.resp, ll=ll, mu=mu, sigma=sigma, 
						pem_pars=pem_pars, a_center_type=a.rater.center, PEM_itermax=PEM_itermax, 
						b.rater.center=b.rater.center, a.rater.center=a.rater.center, 
						a.item.center=a.item.center ) 												
			ll_pem <- res$ll
			pem_parameter_sequence <- res$pem_parameter_sequence
			a.rater <- res$a.rater
			b.rater <- res$b.rater
			a.item <- res$a.item
			tau.item <- res$tau.item
			pi.k <- res$pi.k
			mu <- res$mu
			sigma <- res$sigma				
			PEM <- res$PEM
		}	
		
		#-- save deviance values
		deviance.history[iter+1] <- dev <- -2*ll
		
		
		#-- convergence criteria
		conv <- max( abs(b.rater-b.rater0) , abs( a.rater-a.rater0), abs( tau.item0-tau.item) , abs( a.item - a.item0 ) )
		iter <- iter+1
		devchange <- abs( ( dev - dev0 ) / dev0  )

zz0 <- sirtcat( "  *** calc ll   " , zz0 , active )
	
		#-- print progress			
		res <- rm_facets_print_progress( dev=dev, dev0=dev0, b.rater=b.rater, b.rater0=b.rater0, a.rater=a.rater, a.rater0=a.rater0, 
					tau.item=tau.item, tau.item0=tau.item0, a.item=a.item, a.item0=a.item0, mu=mu, sigma=sigma, 
					iter=iter ) 
	}
				
	#*********
	# arrange OUTPUT
	#---
	ic <- rm_facets_ic( dev=dev, dat2=dat2, VV=VV, RR=RR, maxK=maxK, a.item.center=a.item.center, 
				est.a.item=est.a.item, est.b.rater=est.b.rater, est.a.rater=est.a.rater, 
				est.mean=est.mean, b.rater.center=b.rater.center, a.rater.center=a.rater.center, 
				b.rater.fixed=b.rater.fixed, a.rater.fixed=a.rater.fixed, 
				tau.item.fixed_val=tau.item.fixed_val, a.item.fixed=a.item.fixed ) 	
	#---
	# person
	res <- rm_facets_postproc_person( dat2=dat2, dat2.resp=dat2.resp, procdata=procdata, maxK=maxK, RR=RR, theta.k=theta.k, 
				f.qk.yi=f.qk.yi ) 
	person <- res$person
	EAP.rel <- res$EAP.rel
	
	#---
	# item
	if ( ! is.null(tau.item.fixed) ){
		tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA
		se.tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA		
	}	
    item <- data.frame( "item" = colnames(dat) , 
			"N" = colSums( 1-is.na(dat)) , 
			"M" = colMeans( dat , na.rm=TRUE ) )
	for (kk in 1:K){ 
		item[ , paste0("tau.Cat",kk) ] <- tau.item[,kk] 
	}
    item$a <- a.item

	delta.item <- pcm.conversion(tau.item)$delta
	item$delta <- delta.item
	item$delta_cent <- item$delta - mean( item$delta)

	cat("*********************************\n")
	cat("Item Parameters\n")
	sirt_summary_print_objects(obji=item, digits=3, from=2)	
	
	#---
	# rater
	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
	N <- stats::aggregate( N , list( rater.index ) , sum )[,2]
	M1 <- stats::aggregate( M1 , list( rater.index ) , mean )[,2]	
    rater <- data.frame( "rater" = rater.index1[,1] , "N" = N , "M" = M1 , 	"b" = b.rater ,	"a" = a.rater )
	rater$thresh <- rater$a * rater$b

	#*****
	# dimnames probs
	dimnames(probs)[[1]] <- colnames(dat2)
	#*****
	# expanded item parameters
	ipars.dat2 <- rm_facets_itempar_expanded( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, 
						a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index, 
						theta.k=theta.k, RR=RR ) 
	
	cat("*********************************\n")
	cat("Rater Parameters\n")		
	sirt_summary_print_objects(obji=rater, digits=3, from=2)	
	
	cat("*********************************\n")
	cat("EAP Reliability = " , round(EAP.rel,3) , "\n")		
	
	s2 <- Sys.time()	
    res <-  list( deviance=dev, ic=ic, item=item, rater=rater, person=person, EAP.rel=EAP.rel, mu=mu, 
					sigma=sigma, theta.k=theta.k, pi.k=pi.k, G=1, tau.item=tau.item, se.tau.item=se.tau.item, 
					a.item=a.item, se.a.item=se.a.item, delta.item=delta.item, b.rater=b.rater, 
					se.b.rater=se.b.rater, a.rater=a.rater, se.a.rater=se.a.rater, f.yi.qk=f.yi.qk, 
					f.qk.yi=f.qk.yi, probs=probs, n.ik=n.ik, maxK=maxK, procdata=procdata, iter=iter, s1=s1, s2=s2, 
					tau.item.fixed=tau.item.fixed, item.index=item.index, rater.index=rater.index, 
					ipars.dat2=ipars.dat2, CALL=CALL, deviance.history=deviance.history ) 				
	class(res) <- "rm.facets"
	return(res)
}
