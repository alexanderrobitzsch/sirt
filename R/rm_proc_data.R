## File Name: rm_proc_data.R
## File Version: 0.25

##########################################
# Data preprocessing rater models
rm_proc_data <- function( dat, pid , rater, rater_item_int=FALSE, reference_rater=NULL )
{
	#--- define reference rater
	if ( is.null(reference_rater) ){
	 	reference_rater0 <- sort( paste( rater ))[1]
	}

	#--- item-specific rater parameters and rearrange dataset
	if (rater_item_int){
		if (is.null(reference_rater)){
			reference_rater <- reference_rater0
		}
		res <- rm_proc_create_pseudoraters( dat=dat, rater=rater, pid=pid, reference_rater=reference_rater )
		dat <- res$dat
		rater <- res$rater
		pid <- res$pid
		reference_rater <- res$reference_rater
	}
	
	#-- create rater indices
	rater <- paste( rater )
	# create table of rater indizes
	rater.index <- data.frame( "rater" = sort( unique( rater )) )
	rater.index$rater.id <- seq( 1 , nrow(rater.index) )
	RR <- nrow(rater.index)	
	
	# create table of person indizes
	person.index <- data.frame( "pid" = sort( unique( pid )) )
	person.index$person.id <- seq( 1 , nrow(person.index) )
	PP <- nrow( person.index )
	# number of variables
	VV <- ncol(dat)
	vars <- colnames(dat)    
	# create data frame with crossed items and raters
	dat2 <- data.frame( matrix( NA , nrow=PP , ncol=RR*VV ) )
	colnames(dat2) <- paste0( rep(vars , RR ) , "-" , rep( rater.index$rater , each=VV) )
	rownames(dat2) <- person.index$pid
	
	for (rr in 1:RR){
		ind.rr <- which( rater == rater.index$rater[rr] )
		dat.rr <- dat[ ind.rr  , ]
		pid.rr <- pid[ ind.rr ]
		i1 <- match(  pid.rr , person.index$pid )
		colnames(dat.rr) <- NULL
		dat2[ i1 , VV*(rr-1) + 1:VV ] <- dat.rr
	}
	
	# variable list
	dataproc.vars <- list( "item.index" = rep( 1:VV , RR ), "rater.index" = rep(1:RR , each=VV ) )
	# arrange response data
	dat2.resp <- 1 - is.na(dat2)
	dat20 <- dat2
	dat2[ dat2.resp == 0 ] <- 0
	#--- dat2.ind.resp
	n <- nrow(dat2)
	p <- ncol(dat2)
	K <- max( dat2, na.rm=TRUE) + 1
	dat2.ind.resp <- array( 0 , dim=c(n,p,K) )
	for (kk in 1:K){
		dat2.ind.resp[,,kk] <- dat2.resp * ( dat2 == ( kk - 1 ) )
	}
	
	#--- output
	res <- list( dat2=dat2, dat2.resp=dat2.resp, dat2.NA=dat20, dat=dat, 
				person.index=person.index, rater.index=rater.index, VV=VV, N=PP, RR=RR, 
				dataproc.vars=dataproc.vars, dat2.ind.resp=dat2.ind.resp, rater=rater, pid=pid, dat=dat,
				reference_rater=reference_rater) 
	return(res)
}
#################################################

rm_proc <- rm_proc_data
.prep.data.rm <- rm_proc
