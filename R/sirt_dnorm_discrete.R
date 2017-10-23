## File Name: sirt_dnorm_discrete.R
## File Version: 0.01
## File Last Change: 2017-10-02 15:20:25

sirt_dnorm_discrete <- function(x, mean=0, sd=1)
{
	fx <- stats::dnorm(x, mean=mean, sd=sd)
	fx <- fx / sum(fx)
	return(fx)
}
