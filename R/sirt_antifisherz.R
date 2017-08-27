## File Name: sirt_antifisherz.R
## File Version: 0.01
## File Last Change: 2017-07-11 16:05:55

## see psych::fisherz2r
sirt_antifisherz <- function(z)
{
	( exp(2*z) - 1 ) / ( exp(2*z) + 1 )
}
