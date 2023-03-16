## File Name: lsem_spline.R
## File Version: 0.02

lsem_spline <- function( x, y, method="fmm", n=100)
{
    res <- stats::spline( x=x, y=y, n=n, method=method )
    return(res)
}
