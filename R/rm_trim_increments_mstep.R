## File Name: rm_trim_increments_mstep.R
## File Version: 0.01
## File Last Change: 2017-10-03 15:53:27

rm_trim_increments_mstep <- function( parm, parm0 , max.increment )
{
	increment <- parm - parm0
	increment <- rm_numdiff_trim_increment( increment=increment, max.increment=max.increment, eps2=1E-10) 
	parm <- parm0 + increment
	return(parm)
}
