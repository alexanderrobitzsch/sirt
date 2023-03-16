## File Name: sirt_pem_algorithm_compute_t.R
## File Version: 0.03

sirt_pem_algorithm_compute_t <- function( i, a=1.5, h=0.1)
{
    return( 1 + a^i * h )
}
