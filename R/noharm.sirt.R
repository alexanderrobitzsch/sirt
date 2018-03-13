## File Name: noharm.sirt.R
## File Version: 0.74

########################################
# NOHARM implementation in R
########################################

#---------------
# NOHARM
# Normal-Ogive Harmonic Analysis Robust Method
#---------------
# The other function with linear and loglinear models of parameters
# is denoted as 'noharm.loglm.sirt'.
#---------------

noharm.sirt <- function(dat,weights=NULL,Fval=NULL,Fpatt=NULL,
	Pval=NULL, Ppatt=NULL, Psival=NULL, Psipatt=NULL , dimensions=NULL , 
	lower=rep(0,ncol(dat)), upper=rep(1,ncol(dat)) , wgtm=NULL ,
	modesttype=1 , pos.loading=FALSE , pos.variance = FALSE ,
	pos.residcorr = FALSE , 
	maxiter=1000 , conv=10^(-6) ,  increment.factor=1.01,
	reliability=TRUE ){
	#*****************************************
	# data processing
	e1 <- environment()
	
	ss <- pm0 <- sumwgtm <- model.type <- modtype <- N <- I <- dat0 <- NULL
	wgtm.default <- NULL
	# preprocessing in a subfunction
	res <- .noharm.sirt.preproc( dat , weights , Fpatt , Fval , 
		Ppatt , Pval , Psipatt , Psival, wgtm , dimensions )	
	# attach objects locally
	.attach.environment.sirt( res , envir =e1)
	D <- ncol(Fval)
#	modesttype <- 1		# NOHARM estimation
#	modesttype <- 2		# estimation using tetrachoric correlations
	if (modesttype==2){
		pm2 <- tetrachoric2( dat , method="Bo")
		pm <- pm2$rho
		betaj <- pm2$tau
		lower <- 0*lower
		upper <- 1+lower
	}
	#----
	# compute betaj
	if (modesttype==1){
		betaj <- - stats::qnorm( ( diag(pm) - lower ) / ( upper - lower ) )
	}
	# compute bjk coefficients
	# include lower and upper asymptotes here
	b0 <- lower + (upper-lower) * stats::pnorm( - betaj )
	b1 <- (upper-lower) * stats::dnorm( betaj )
	b2 <- (upper-lower) * betaj * stats::dnorm( betaj ) / sqrt(2)
	b3 <- (upper-lower) * ( betaj^2 -1 ) * stats::dnorm(betaj) / sqrt(6)
	# create fixed cofficients
	b0.jk <- as.matrix( TAM::tam_outer(b0,b0) )
	b1.jk <- as.matrix( TAM::tam_outer(b1,b1) )
	b2.jk <- as.matrix( TAM::tam_outer(b2,b2) )
	b3.jk <- as.matrix( TAM::tam_outer(b3,b3) )

	parchange <- 1
	changeF <- changeP <- changePsi <- 0
	iter <- 0
	maxincrement <- .2
	estPsi <- 1 * ( sum( Psipatt > 0 ) > 0 )
	estF <- 1 * ( sum( Fpatt > 0 ) > 0 )	
	estP <- 1 * ( sum( Ppatt > 0 ) > 0 )	
	estpars <- list("estF"=estF, "estP"=estP,"estPsi"=estPsi)
	eps <- 2*conv
	
	#*******************
	# algorithm
	while( ( iter < maxiter ) & ( parchange > conv ) ){
		maxincrement <- maxincrement  / increment.factor
		#---- update F
		if (estF==1){	
			res <- noharm_estFcpp( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I, D=D, 
						b0jk=b0.jk, b1jk=b1.jk, b2jk=b2.jk, b3jk=b3.jk, wgtm=wgtm, pm=pm, Psival=Psival, 
						Psipatt=Psipatt, maxincrement=maxincrement, modtype=modesttype ) 
			changeF <- res$change	
			Fval <- res$Fval_	
			if ( pos.loading ){	Fval[ Fval < 0 ] <- eps  }
		}
		#---- update P
		if (estP==1){
			Pval_old00 <- Pval
			res <- noharm_estPcpp( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I, D=D, 
						b0jk=b0.jk, b1jk=b1.jk, b2jk=b2.jk, b3jk=b3.jk, wgtm=wgtm, pm=pm, Psival=Psival, 
						Psipatt=Psipatt, maxincrement=maxincrement, modtype=modesttype )
			changeP <- res$change
			Pval <- res$Pval_		
			if ( pos.variance ){	
				diag(Pval)[ diag(Pval) < 0 ] <- eps	
			}			
			diag(Pval)[ is.na(diag(Pval)) ] <- eps
		}
		#---- update Psi
		if (estPsi==1){
			res <- noharm_estPsicpp( Fval=Fval, Pval=Pval, Fpatt=Fpatt, Ppatt=Ppatt, I=I, D=D, 
						b0jk=b0.jk, b1jk=b1.jk, b2jk=b2.jk, b3jk=b3.jk, wgtm=wgtm, pm=pm, Psival=Psival, 
						Psipatt=Psipatt, maxincrement=maxincrement, modtype=modesttype ) 
			changePsi <- res$change
			Psival <- res$Psival_	
			if ( pos.residcorr ){	Psival[ Psival < 0 ] <- eps  }
		}
		parchange <- max( c(changeP,changeF,changePsi) )
		iter <- iter + 1 
	}
	#--------------------------------

	#****************
	# calculate final constants
	if (modesttype==2){
		dj2 <- diag( Fval %*% Pval %*% t(Fval) ) 
		Fval <- Fval / sqrt( 1 - dj2 )
	}
	# recalculation of f0 coefficient (final constant)
	dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )
	ej <- sqrt( 1 + dj^2 )
	f0 <- - betaj * ej	
	# uniquenesses
	uqn <- 1 - ( dj^2 / ( 1 + dj^2 ) )
	# standardized loadings
	loadingsF <- Fval / ej
	#****************
	# post processing
	residuals <- .noharm.est.residuals( Fval , Pval , Fpatt , Ppatt , 
		I , D ,  b0.jk , b1.jk , b2.jk , b3.jk , wgtm , pm ,
		Psival , Psipatt )	
	rmsr <- sqrt( sum( residuals^2 * wgtm ) / sum( wgtm) )
	RMW <- residuals * wgtm 
	PMW <- pm * wgtm
	tanaka <- 1 - sum( diag( RMW %*% RMW ) ) /  sum(diag( PMW %*% PMW ))	
	# some labels
	rownames(Psival) <- colnames(Psival) <- colnames(dat)
	if ( modesttype==2){ tanaka <- rmsr <- NA }
	
	#****************
	# output
	res <- list( "tanaka"=tanaka , "rmsr"=rmsr , "N.itempair"=ss , 
				"pm"=pm0 , "wgtm"=wgtm , "sumwgtm" = sumwgtm , 
				"lower"=lower , "upper"=upper ,
				"residuals" = residuals , "final.constants"=f0 ,
				"factor.cor"=Pval , "thresholds" = betaj , 
				"uniquenesses"=uqn , "loadings"=loadingsF , 
				"loadings.theta"=Fval , "residcorr" = Psival , 
				"model.type" = model.type , "modtype"=modtype , "Nobs"=N , 
				"Nitems" = I , 
				"Fpatt"=Fpatt , "Ppatt"=Ppatt, "Psipatt"=Psipatt ,
				"dat"=dat0 , 	
				"systime"=Sys.time() , "dimensions"=D , 
				"display.fit"=5,
				"iter"=iter )
	
	#*********
	# number of estimated parameters
	Nestpars <- list("total"=0)
	Nestpars$total <- Nestpars$thresh <- I
	Nestpars$F <- sum( Fpatt > 0 )
	Nestpars$total <- Nestpars$total + Nestpars$F
	Nestpars$P <- sum( diag(Ppatt) == 1 ) / 2 +  sum( Ppatt == 1 ) / 2 
	Nestpars$total <- Nestpars$total + Nestpars$P
	Nestpars$Psi <- 1/2 * sum( Psipatt == 1 )
	Nestpars$total <- Nestpars$total + Nestpars$Psi
	res$Nestpars <- Nestpars
	#******************
    # chi square statistics	
	if (modesttype==1){	
		RM <- residuals
		PV <- diag( pm ) 
		g1 <- sqrt( outer( PV * (1-PV) , PV * ( 1-PV ) ) )
		rM <- ( RM / g1  )
		zM <- 0.5 * log( 1 + rM ) - 0.5 * log( 1 - rM ) 
		# chi square
		res$chisquare <- X2 <- ( N - 3 ) * sum( zM^2 )		
		# res$df <- df <- 0.5*I*(I+1) - Nestpars$total	
		res$df <- df <- I + sumwgtm - Nestpars$total
		res$chisquare_df  <- res$chisquare / res$df
		# calculate RMSEA
		res$rmsea <- rmsea <- sqrt(max(c( (X2 / res$Nobs ) / df - 1/res$Nobs , 0)))
		# calculate p values
		res$p.chisquare <- 1 - pchisq( res$chisquare , df = res$df )
	}

	#************************
	# Green-Yang reliability
	v1 <- 1 * ( ( wgtm - diag(wgtm) ) > 0 )
	# positive definiteness
    res$omega.rel <- NA
		L0 <- sqrt( diag(Pval ) )
		N1 <- length(L0)
		L1inv <- L1 <- matrix( 0 , N1 , N1 )
		rownames(L1) <- rownames(L1inv) <- colnames(L1) <- 
				colnames(L1inv) <- rownames(Pval)
		diag(L1) <- L0	
		diag(L1inv) <- 1 / L0
		Fval2 <- Fval %*% L1
		Pval2 <- L1inv %*% Pval %*% L1inv
		dj <- sqrt( diag( Fval2 %*% Pval2 %*% t(Fval2) ) )
		ej <- sqrt( 1+dj^2 )
		Fval2 <- Fval2 / ej
		standardized.solution <- list( "Fval"=Fval2 , "Pval" = Pval2 )					
	if (reliability){
		res$omega.rel <- reliability.nonlinearSEM(facloadings= Fval2 , 
			thresh=res$thresholds, resid.cov = res$residcorr  , cor.factors = Pval2 )$omega.rel
	}
	if ( sum(v1) + I < I^2 ){
		res$omega.rel <- NA  
	} 
	
	#***********************
	# rotated solution	
	if (model.type=="EFA"){	
		m1 <- stats::promax(res$loadings)		
		p1 <- matrix( 0 , nrow=I,D)
		for (dd in 1:D){ p1[,dd] <- m1$loadings[,dd] }
		colnames(p1) <- colnames(Pval)
		rownames(p1) <- colnames(dat)
		res$promax <- p1
		res$factor.cor <- solve( crossprod( m1$rotmat ) ) 
		rownames(res$factor.cor) <- colnames(res$factor.cor) <- colnames(Pval)
		# conversion to THETA parametrization
		h1 <- rowSums( p1^2 )
		p2 <- p1 / sqrt( max(1 - h1 , .0001 ) )	
		res$promax.theta <- p2
				}
    if (modesttype==2){ res$tetracor <- pm }
				
	#***----- more output
	res$estpars <- estpars	
	res$modesttype <- modesttype
	res$guesses <- res$lower
	res$wgtm.default <- wgtm.default
	res$standardized.solution <- standardized.solution
	#****************************
	# return
	class(res) <- "noharm.sirt"
	return(res)
}
##################################################


