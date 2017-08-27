## File Name: mirt.wrapper.posterior.R
## File Version: 0.24
## File Last Change: 2017-03-17 18:01:23

#############################################################
# calculation of posterior
mirt.wrapper.posterior <- function( mirt.obj , weights=NULL ){
	# adapt function for multiple groups
    # ****
	mobj <- mirt.obj
	# extract theta
	# Theta <- mobj@Theta
	Theta <- mobj@Model$Theta
	# extract theta distribution
	# pi.k <- mobj@bfactor$Prior[[1]]
	# pi.k <- mobj@Prior[[1]]
	pi.k <- mobj@Internals$Prior[[1]]
	# extract full data
	fulldata <- mobj@Data$fulldata[[1]]
	#load items
	I <- ncol( mobj@Data$data )
	items <- vector('list', I)
	for(ii in 1:I){
		items[[ii]] <- mirt::extract.item(mobj, ii)
	}				
	# check whether prodlist exists
    # prodlist <- attr(mobj@pars, "prodlist") 
	prodlist <- mobj@Model$prodlist
	# Theta <- mobj@Theta 
    Theta1 <- Theta
	if ( length(prodlist) > 0 ){
		Theta1 <- mirt_prodterms(Theta, prodlist)	
	}				
	# item-wise probabilities for each Theta
	traces <- Probtrace_sirt(items, Theta1)	
	# log-Likelihood
	f.yi.qk <- exp( fulldata %*% t(log(traces))  )
	# compute individual posterior
	N <- nrow( fulldata )
	TP <- length(pi.k)
	piM <- matrix( pi.k , nrow=N , ncol=TP , byrow=TRUE )
	f.qk.yi <- f.yi.qk * piM 
	f.qk.yi <- f.qk.yi / matrix( rowSums( f.qk.yi ) , nrow=N , ncol=TP , byrow=FALSE )
	# maximum category
	maxK <- apply( mobj@Data$data , 2 , max , na.rm=TRUE)+1
	resp.ind <- 1 - is.na(mobj@Data$data)
	resp <- mobj@Data$data
	resp[ resp.ind == 0 ] <- 0
	# calc counts
	group <- NULL	# only applies to single groups for now
	if (is.null(weights) ){ 
		pweights <- rep(1,N) 
	} else {
		pweights <- weights 
	}
    # Theta is only used for calculating dimension size					
	n.ik <- mirt.wrapper.calc.counts( resp, theta=Theta , resp.ind=resp.ind , 
				group=group , maxK=max(maxK) , pweights=pweights , hwt=f.qk.yi )
    probs <- traces		
    probs <- array( probs , dim = c(TP,max(maxK),I) )	
	probs <- aperm( probs , c(3,2,1) )		
	# result list
	res <- list( "theta.k" = Theta , "pi.k" = pi.k ,
			"f.yi.qk" = f.yi.qk , "f.qk.yi" = f.qk.yi ,
			"n.ik"=n.ik , "probs" = probs ,
			"N"= N , "TP"=TP , "I" = I , "data" = mobj@Data$data ,
			"maxK" = maxK )
	class(res) <- "mirt"
	return(res)
}

