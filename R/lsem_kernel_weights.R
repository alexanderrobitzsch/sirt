## File Name: lsem_kernel_weights.R
## File Version: 0.05

lsem_kernel_weights <- function(x, x0, bw, kernel="gaussian")
{
    if (kernel=="gaussian"){
        wgt <- exp( - (x - x0)^2 / (2*bw^2) )
    }
    if (kernel=="uniform"){
        wgt <- 1*(abs(x-x0) <=bw)
    }
    if (kernel=="epanechnikov"){
        z <- (x-x0)/bw
        wgt <- abs(z<1)*3/4*(1-z^2)
    }
    return(wgt)
}
