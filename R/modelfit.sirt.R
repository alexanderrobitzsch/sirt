## File Name: modelfit.sirt.R
## File Version: 1.02


# model fit in sirt
modelfit.sirt <- function( object )
{
	#****
	# object of class tam.mml, tam.mml.2pl or tam.fa
	# Note that only dichotomous responses are allowed
	if (class(object) %in% c("tam.mml","tam.mml.2pl")){
		mod <- object
		posterior <- mod$hwt
		probs <- mod$rprobs
		dat <- mod$resp
		dat[ mod$resp.ind == 0 ] <- NA
	}
	#*****
	# rasch.mml
	if (class(object)=="rasch.mml"){
		mod <- object
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
		probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
		probs[ , 2 , ] <- t(prob1)
		probs[ , 1 , ] <- 1 - t(prob1)
		dat <- mod$dat
	}
	#*****
	# rasch.mirtlc
	if (class(object)=="rasch.mirtlc"){
		mod <- object$estep.res
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
		probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
		probs[ , 2 , ] <- t(prob1)
		probs[ , 1 , ] <- 1 - t(prob1)
		dat <- object$dat
	}	
	#******
	# rasch.pml
	if ( class(object) !="rasch.pml"){ pmlobject <- NULL } else {
		data <- NULL ; posterior <- NULL ; probs <- NULL ; pmlobject <- object }	
	#*******
	# smirt	
	if (class(object) == "smirt"){
		# note that for polytomous response data some adaptations are
		# necessary: see modelfit in the CDM package
		mod <- object
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
	}
	#*******
	# smirt	
	if (class(object) == "gom"){
		mod <- object
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
	}					
	#*******
	# rm.facets
#	if (class(object) %in% c("rm.facets") ){
#		mod <- object
#		probs <- mod$probs
#		posterior <- mod$f.qk.yi
#		dat <- mod$procdata$dat2.NA
#	}						
					
	#*******
	# mirt	
	if (class(object) == "ConfirmatoryClass" | class(object)=="ExploratoryClass" ){
		mod <- object
		mod <- mirt.wrapper.posterior(mod)		
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
	}					
	#******
	# R2noharm, noharm.sirt
    if ( class(object) %in% c("R2noharm","noharm.sirt") ){
		# exclusion criteria for noharm.sirt
		if ( class(object) == "noharm.sirt") {
			if ( object$estpars$estPsi > 0 ){
				stop("Model fit cannot be calculated because of correlated residuals")
			}
			if ( ! ( object$wgtm.default ) ){
				stop("Model fit cannot be calculated because not all item pairs are used for estimation")
			}										
		}
		# evaluation of posterior
		mod <- R2noharm.EAP(noharmobj=object, theta.k = seq(-6, 6, len = 15 ) ,
				print.output=FALSE )
		probs <- aperm( mod$probs , c(1,3,2) )
		posterior <- mod$posterior
		dat <- object$dat		  		  
	}
	# calculate modelfit.cor
    if ( class(object) == "rasch.pml" ){
		res <- modelfit.cor.sirt.pml( data = dat , posterior =posterior , probs = probs ,
				           pmlobject=pmlobject)
	} else {
		res <- CDM::modelfit.cor2( data = dat , posterior =posterior , probs = probs )								
	}
    class(res) <- "modelfit.sirt"							
	return(res)
}
################################################################################
