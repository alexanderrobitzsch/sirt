## File Name: rm_numdiff_discrete_differences.R
## File Version: 0.03

rm_numdiff_discrete_differences <- function(ll0, ll1, ll2, h)
{
	# first derivative
	# f(x+h)-f(x-h) = 2*f'(x)*h
	d1 <- ( ll1 - ll2  ) / ( 2 * h )
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
	d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
	#--- output
	res <- list(d1=d1, d2=d2)
	return(res)
}
