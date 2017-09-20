## File Name: sirt_round_vector.R
## File Version: 0.01
## File Last Change: 2017-09-20 10:26:56


sirt_round_vector <- function(x , digits)
{
	if (is.numeric(x)){
		x <- round( x , digits )
	}
	return(x)
}
