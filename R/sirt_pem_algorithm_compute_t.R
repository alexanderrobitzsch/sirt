## File Name: sirt_pem_algorithm_compute_t.R
## File Version: 0.01
## File Last Change: 2017-10-03 12:23:47

sirt_pem_algorithm_compute_t <- function( i, a=1.5, h=0.1)
{
	return( 1 + a^i * h )
}
