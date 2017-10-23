## File Name: rm_squeeze.R
## File Version: 0.01
## File Last Change: 2017-10-02 14:50:56

rm_squeeze <- function(x, lower, upper )
{
	x[ x < lower ] <- lower
	x[ x > upper ] <- upper
	return(x)
}
