## File Name: lsem_kernel_weights.R
## File Version: 0.01

lsem_kernel_weights <- function(x, x0, bw)
{
    wgt <- exp( - (x - x0)^2 / (2*bw^2) )
    return(wgt)
}
