## File Name: rm_proc_create_pseudoraters.R
## File Version: 0.02

rm_proc_create_pseudoraters <- function( dat, rater, pid, reference_rater=NULL )
{
	dat0 <- dat
	rater0 <- rater
	pid0 <- pid
	items <- colnames(dat)
	dat <- NULL
	pid <- NULL
	rater <- NULL
	I <- length(items)
	for (ii in 1:I){
		dat_ii <- NA*dat0
		dat_ii[, ii ] <- dat0[,ii]
		rater_ii <- sort(paste0( rater0 , "-" , items[ii] ))
		pid_ii <- pid0
		rater <- c( rater, rater_ii )
		dat <- rbind( dat, dat_ii)
		pid <- c( pid, pid_ii)
	}
	if ( ! is.null(reference_rater) ){
		reference_rater <- paste0(reference_rater, "-" , items )
	}	
	#-- output
	res <- list( dat=dat, rater=rater, pid=pid, reference_rater=reference_rater)
	return(res)
}
