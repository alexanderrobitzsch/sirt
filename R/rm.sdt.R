## File Name: rm.sdt.R
## File Version: 8.553

#################################################################
# Hierarchical rater model
# MML estimation
rm.sdt <- function( dat , pid , rater ,Qmatrix=NULL , theta.k=seq(-9,9,len=30) , 	
	est.a.item=FALSE , est.c.rater= "n" , 
	est.d.rater= "n" , est.mean = FALSE , skillspace="normal" , 
	tau.item.fixed = NULL , a.item.fixed = NULL , 
	d.min=.5 , d.max=100 ,  d.start = 3,  
	d.prior = c(3,100), c.prior=c(3,100), tau.prior=c(0,1000), a.prior=c(1,100), 
	max.increment=1 , numdiff.parm=.00001 , maxdevchange=.10 ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001, fac_incr=.99,
	PEM=FALSE, PEM_itermax=maxiter )
{
	#..........................................................
	CALL <- match.call()	
	s1 <- Sys.time()
	theta.k0 <- theta.k
	pi.k <- sirt_dnorm_discrete(x=theta.k, mean=0, sd=1)
	max.b.increment <- max.increment		
	a_center_type <- 2	
	
	#-- process data
	procdata <- res <- rm_proc_data( dat=dat , rater=rater , pid=pid )
	dat2 <- as.matrix(res$dat2)
	dat2.resp <- as.matrix(res$dat2.resp)
	rater.index1 <- res$rater.index
	dataproc.vars <- res$dataproc.vars
	VV <- res$VV
	RR <- res$RR
	item.index <- res$dataproc.vars$item.index 
	rater.index <- res$dataproc.vars$rater.index 
	dat2.ind.resp <- res$dat2.ind.resp	
	
	deviance.history <- rep(NA, maxiter )
	
	# maximum categories
	maxK <- sirt_colMaxs(dat)
	K <- max( maxK )
	if ( is.null(Qmatrix) ){
		Qmatrix <- matrix( 1:K , nrow=VV , ncol=K , byrow=TRUE)
	}
	TP <- length(theta.k)
	I <- VV*RR
	
	# define constraints on tau.item parameters
	# if not all categories are observed
	if ( is.null( tau.item.fixed )){	
		tau.item.fixed <-  rm_determine_fixed_tau_parameters( K=K, maxK=maxK, VV=VV ) 
	}
	
	# starting values for item difficulties
	tau.item <- matrix( 0 , nrow=VV , ncol=K )
	rownames(tau.item) <- colnames(dat)
	
	tau.item <- matrix( seq( -2 , 2 , len=K ) , nrow=VV , ncol=K , byrow=TRUE )
	if ( ! is.null(tau.item.fixed) ){
	    tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
	}
		
	a.item <- rep( 1 , VV )
	if (skillspace == "discrete" ){
		est.mean <- TRUE
	}
	if ( ! is.null( a.item.fixed ) ){
		est.a.item <- TRUE
		a.item[ a.item.fixed[,1] ] <- a.item.fixed[,2]
	}
	
	# rater parameter
	d.rater <- matrix( d.start , nrow=I , ncol=1 )
	c.rater <- matrix( d.start*((1:K) - .5 ) , nrow=I , ncol=K , byrow=TRUE )
	
	# set c.rater for fixed items to 99
	c.rater.fixed <- NULL
	if ( ! is.null( tau.item.fixed ) ){
		tau1 <- tau.item.fixed[ tau.item.fixed[,3] == 99 , , drop=FALSE]		
		ind <- match( item.index , tau1[,1] )
		c.rater.fixed <- tau1[ ind , ]
        c.rater.fixed[,1] <- seq( 1 , nrow(c.rater.fixed) )	
		c.rater.fixed[,3] <- 999
		c.rater.fixed <- c.rater.fixed[ ! is.na( c.rater.fixed[,2] ) , ] 
		c.rater[ c.rater.fixed[,1:2] ] <- c.rater.fixed[,3]		
	}
	
	#--- indices for derivatives	
    diffindex <- rm_sdt_prepare_diffindex( item.index=item.index, rater.index=rater.index, I=I, est.c.rater=est.c.rater, 
					est.d.rater=est.d.rater ) 
		
    # init standard errors
	se.d.rater <- NA*d.rater
	se.c.rater <- NA*c.rater
	se.a.item <- NA*a.item
	
	d.rater.incr <- 2 
	tau.item.incr  <- max.b.increment
	c.rater.incr <- max.b.increment
	a.item.incr <- max.b.increment
	
	#-- preliminaries PEM acceleration
	if (PEM){
		res <- rm_sdt_pem_inits( tau.item=tau.item, a.item=a.item, c.rater=c.rater, d.rater=d.rater, 
					est.a.item=est.a.item, a.item.fixed=a.item.fixed, skillspace=skillspace, PEM=PEM ) 
		PEM <- res$PEM
		pem_pars <- res$pem_pars
		center_log_a <- res$center_log_a
		pem_parameter_index <- res$pem_parameter_index
		pem_parameter_sequence <- res$pem_parameter_sequence
	}
	
	# inits
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	mu <- 0
	sigma <- 1
	disp <- "...........................................................\n"	
	prob.item <- NULL
	prob.rater <- NULL
	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter ) ){
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
		
		# previous values
		d.rater0 <- d.rater
		tau.item0 <- tau.item
		dev0 <- dev
		mu0 <- mu
		sigma0 <- sigma
		a.item0 <- a.item
		c.rater0 <- c.rater

		# calculate probabilities
		res <- rm_hrm_calcprobs( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
					d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR ) 
        probs <- res$prob.total				
		prob.rater <- res$prob.rater
		prob.item <- res$prob.item
		
		#-- calculate posterior
		res <- rm_posterior( dat2=dat2, dat2.resp=dat2.resp, TP=TP, pi.k=pi.k, K=K, I=I, probs=probs, dat2.ind.resp=dat2.ind.resp ) 
		f.yi.qk <- res$f.yi.qk
		f.qk.yi <- res$f.qk.yi
		n.ik <- res$n.ik
		N.ik <- res$N.ik
		pi.k <- res$pi.k
		ll <- res$ll
																					
    	#-- estimate tau.item parameters
		res <- rm_hrm_est_tau_item( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
					d.rater=d.rater, item.index=item.index, rater.index=rater.index, n.ik=n.ik, 
					numdiff.parm=numdiff.parm, max.b.increment=tau.item.incr, theta.k=theta.k, 
					msteps=msteps, mstepconv=mstepconv, tau.item.fixed=tau.item.fixed, 
					prob.rater=prob.rater, tau.item0=tau.item0, tau.prior=tau.prior ) 
		tau.item <- res$tau.item
		se.tau.item <- res$se.tau.item
		g1  <- abs( tau.item0 - tau.item )
		# tau.item.incr <- ifelse( tau.item.incr > g1 , g1 , tau.item.incr )						
		tau.item.incr <- fac_incr^iter
		prob.item <- res$prob.item
		
		#-- estimate a.item parameter
		if (est.a.item){
			res <- rm_hrm_est_a_item( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
						d.rater=d.rater, item.index=item.index, rater.index=rater.index, n.ik=n.ik, 
						numdiff.parm=numdiff.parm, max.b.increment=a.item.incr, theta.k=theta.k, msteps=msteps, 
						mstepconv=mstepconv, prob.rater=prob.rater, a.item.fixed=a.item.fixed, a_center_type = a_center_type,
						a.item0=a.item0, a.prior=a.prior) 
			a.item <- res$a.item
			se.a.item <- res$se.a.item
			# a.item.incr <- abs( a.item0 - a.item )
			a.item.incr <- fac_incr^iter
			prob.item <- res$prob.item
		}
				
		#-- estimate d.rater parameter
		if (est.d.rater!="n"){
			res <- rm_hrm_est_d_rater( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
						d.rater=d.rater, item.index=item.index, rater.index=rater.index, n.ik=n.ik, 
						numdiff.parm=numdiff.parm, max.b.increment=d.rater.incr, theta.k=theta.k, 
						msteps=msteps, mstepconv=mstepconv, d.min=d.min, d.max=d.max, est.d.rater=est.d.rater, 
						prob.item=prob.item, d.rater0=d.rater0, diffindex = diffindex$d.rater, d.prior=d.prior ) 
			d.rater <- res$d.rater
			se.d.rater <- res$se.d.rater
			g1 <- fac_incr^iter
			d.rater.incr <- max(g1)			
		}				
					
		#-- estimate c.rater parameter
		if( est.c.rater != "n" ){			
			res <- rm_hrm_est_c_rater( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
						d.rater=d.rater, item.index=item.index, rater.index=rater.index, n.ik=n.ik, 
						numdiff.parm=numdiff.parm, max.b.increment=c.rater.incr, theta.k=theta.k, 
						msteps=msteps, mstepconv=mstepconv, est.c.rater=est.c.rater, prob.item=prob.item, 
						c.rater.fixed=c.rater.fixed, c.rater0=c.rater0, diffindex=diffindex$c.rater, c.prior=c.prior ) 
			c.rater <- res$c.rater
			se.c.rater <- res$se.c.rater
			g1 <- abs( c.rater0 - c.rater )
			# g1 <- .99^iter*abs( c.rater0 - c.rater )
			g1 <- fac_incr^iter
			c.rater.incr <- max(g1)
		}						
			
		#-- update distribution
		res <- rm_smooth_distribution( theta.k=theta.k, pi.k=pi.k, est.mean=est.mean, skillspace=skillspace ) 
		pi.k <- res$pi.k
		mu <- res$mu
		sigma <- res$sigma
				
		#-- PEM acceleration
		if (PEM){
			res <- rm_sdt_pem_acceleration( iter=iter, pem_parameter_index=pem_parameter_index, 
						pem_parameter_sequence=pem_parameter_sequence, c.rater=c.rater, Qmatrix=Qmatrix, 
						tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, d.rater=d.rater, 
						item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR, dat2=dat2, 
						dat2.resp=dat2.resp, pi.k=pi.k, dat2.ind.resp=dat2.ind.resp, center_log_a=center_log_a, 
						ll=ll, mu=mu, sigma=sigma, pem_pars=pem_pars, a_center_type=a_center_type,
						PEM_itermax=PEM_itermax ) 
			ll_pem <- res$ll
			pem_parameter_sequence <- res$pem_parameter_sequence
			c.rater <- res$c.rater
			d.rater <- res$d.rater
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
		conv <- max( abs(c.rater-c.rater0) , abs( c.rater-c.rater0) , 
					abs( tau.item0-tau.item) , abs( a.item - a.item0 ) )
		iter <- iter+1
		devchange <- abs( ( dev - dev0 ) / dev0  )

		#-- print progress			
		res <- rm_sdt_print_progress( dev=dev, dev0=dev0, c.rater=c.rater, c.rater0=c.rater0, d.rater=d.rater, d.rater0=d.rater0, 
					tau.item=tau.item, tau.item0=tau.item0, a.item=a.item, a.item0=a.item0, mu=mu, sigma=sigma, iter=iter ) 
		
	}
	#---------------------------- end EM algorithm

	# *********
	# arrange OUTPUT
	
	# c parameters
	if ( ! is.null( c.rater.fixed ) ){
		c.rater[ c.rater.fixed[,1:2] ] <- NA
	}
	
	#---
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$VV <- VV
	ic$RR <- RR
	#****
	# skill parameters
	ic$np.skill <- 0
	if ( skillspace == "normal" ){
		ic$np.skill <- 1 + est.mean
	}
	if ( skillspace == "discrete" ){
		ic$np.skill <- length(theta.k) - 1
	}	
	#*****
	# item parameters
	ic$np.item <- VV * max( maxK )
	if ( ! is.null(tau.item.fixed)){
		ic$np.item <- ic$np.item - nrow(tau.item.fixed)
	}
	if ( is.null( a.item.fixed ) ){ 
		ic$np.item <- ic$np.item + est.a.item*(VV-1)
	} else {
		ic$np.item <- ic$np.item + VV - 1
	}
						
	#*****
	# rater parameters
	ic$np.rater <- 0
	if ( est.d.rater=="e" ){ ic$np.rater <- ic$np.rater + 1 }
	if ( est.d.rater=="i" ){ ic$np.rater <- ic$np.rater + VV }	
	if ( est.d.rater=="r" ){ ic$np.rater <- ic$np.rater + RR }	
	if ( est.d.rater=="a" ){ ic$np.rater <- ic$np.rater + I }		
	if ( est.c.rater=="e" ){ 
		ic$np.rater <- ic$np.rater + K 
	}
	if ( est.c.rater=="i" ){ 
		ic$np.rater <- ic$np.rater + sum(maxK)
	}	
	if ( est.c.rater=="r" ){ 
		ic$np.rater <- ic$np.rater + K*RR 
	}		
	if ( est.c.rater=="a" ){ 
		ic$np.rater <- ic$np.rater + K*I 
		if ( ! is.null( c.rater.fixed ) ){
			ic$np.rater <- ic$np.rater - nrow( c.rater.fixed )
		}			
	}		
	ic$np <- ic$np.skill + ic$np.item + ic$np.rater
    
	#-- compute information criteria
	ic <- rm_ic_criteria(ic=ic)	
	
	#---
	# person parameters
	res <- rm_facets_postproc_person( dat2=dat2, dat2.resp=dat2.resp, procdata=procdata, maxK=maxK, RR=RR, theta.k=theta.k, 
				f.qk.yi=f.qk.yi ) 
	person <- res$person
	EAP.rel <- res$EAP.rel					
	
	#---
	# item
	if (!is.null(tau.item.fixed)){
		K <- max(maxK)
        I <- nrow(tau.item)
		for (ii in 1:I){
			if ( maxK[ii] < K ){
				for (kk in seq(maxK[ii]+1,K) ){
					tau.item[ii,kk] <- NA
				}
			}
		}
		se.tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA								
	}
							
    item <- data.frame( "item" = colnames(dat) , 
			"N" = colSums( 1-is.na(dat)) , 
			"M" = colMeans( dat , na.rm=TRUE ) )
	for (kk in 1:K){ item[ , paste0("tau.Cat",kk) ] <- tau.item[,kk] }
    item$a <- a.item
#	item$b <- rowMeans(tau.item)

	# latent mean and standard deviation
	me1 <- rep(NA,VV)
	sd1 <- rep(NA,VV)
	for (ii in 1:VV){
		pii <- prob.item[ii,,]
		qii <- matrix( c(0,Qmatrix[ii,]) , nrow= K+1 , ncol=ncol(pii) )
		me1[ii] <- sum( colSums( qii * pii ) * pi.k )
		sd1[ii] <- sqrt( sum( colSums( qii^2 * pii ) * pi.k ) - me1[ii]^2  )
	}
	item$latM <- me1
	item$latSD <- sd1	
	
	
	cat("*********************************\n")
	cat("Item Parameters\n")
	sirt_summary_print_objects(obji=item, digits=3, from=2)		
	
	#---
	# rater
	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
    rater <- data.frame( "item.rater" = colnames(dat2) , 
			"N" = N , "M" = M1 , "d" = d.rater )
    for (zz in 1:(ncol(c.rater) ) ){
		rater[ , paste0("c_",zz)] <- c.rater[,zz] 
	}
	# transformed c parameters
    for (zz in 1:(ncol(c.rater) ) ){
		# rater[ , paste0("c_",zz,".trans")] <- c.rater[,zz] / d.rater / K 
        rater[ , paste0("c_",zz,".trans")] <- c.rater[,zz] / d.rater	
	}

	rater <- rater[ order( paste( rater$item.rater) ) , ]
	rownames(rater) <- NULL
	rownames(item) <- NULL

	rt1 <- rater
	l1 <- paste(rt1$item.rater)
	l2 <- strsplit( l1 , split="-" )
	rt1$item <- unlist( lapply( l2 , FUN = function(uu){ uu[[1]] } ) )
	rt1$rater <- unlist( lapply( l2 , FUN = function(uu){ uu[[2]] } ) )
	
	#*****
	# dimnames probs
	dimnames(probs)[[1]] <- colnames(dat2)
	#*****			
						
	#*****
	# distribution	
	skill.distribution <- data.frame("theta.k" = theta.k , "pi.k" = pi.k )
	
	#*****
	# labels
	dimnames(prob.item) <- list( colnames(dat) , paste0("Cat" , 0:K) , NULL )
	dimnames(prob.rater) <- list( colnames(dat2) , paste0("Cat" , 0:K) , NULL )	
						
	cat("*********************************\n")
	cat("Rater Parameters\n")		
	sirt_summary_print_objects(obji=rater, digits=3, from=2)	
	
	cat("*********************************\n")
	cat("EAP Reliability = " , round(EAP.rel,3) , "\n")		
	
	s2 <- Sys.time()

    res <- list( deviance=dev, ic=ic, item=item, rater=rater, person=person, EAP.rel=EAP.rel, mu=mu, 
					sigma=sigma, theta.k=theta.k, pi.k=pi.k, G=1, tau.item=tau.item, se.tau.item=se.tau.item, 
					a.item=a.item, se.a.item=se.a.item, c.rater=c.rater, se.c.rater=se.c.rater, 
					d.rater=d.rater, se.d.rater=se.d.rater, f.yi.qk=f.yi.qk, f.qk.yi=f.qk.yi, probs=probs, 
					prob.item=prob.item, n.ik=n.ik, maxK=maxK, pi.k=pi.k, procdata=procdata, iter=iter, 
					theta.k=theta.k, Qmatrix=Qmatrix, s1=s1, s2=s2, tau.item.fixed=tau.item.fixed, rater2=rt1, 
					maxK=maxK, skill.distribution=skill.distribution, skillspace=skillspace, CALL=CALL,
					deviance.history=deviance.history) 
	class(res) <- "rm.sdt"
	return(res)
}

