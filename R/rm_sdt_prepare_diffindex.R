## File Name: rm_sdt_prepare_diffindex.R
## File Version: 0.02

rm_sdt_prepare_diffindex <- function( item.index, rater.index, I, est.c.rater, est.d.rater )
{
	res <- NULL
	
	#---- est.c.rater
	if (est.c.rater=="r"){ diffindex <- rater.index }
	if (est.c.rater=="i"){ diffindex <- item.index }
	if (est.c.rater=="e"){ diffindex <- rep(1,I) }	
	if (est.c.rater=="a"){ diffindex <- 1:I }		
	res$c.rater <- diffindex
	
	#---- est.d.rater
	if (est.d.rater=="r"){ diffindex <- rater.index }
	if (est.d.rater=="i"){ diffindex <- item.index }
	if (est.d.rater=="e"){ diffindex <- rep(1,I) }	
	if (est.d.rater=="a"){ diffindex <- 1:I }		
	res$d.rater <- diffindex
	
	return(res)
}
