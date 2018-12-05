## File Name: sirt_eigenvalues.R
## File Version: 0.21


# calculation of first D eigenvalues
sirt_eigenvalues <- function(X, D, maxit=200, conv=10^(-6) )
{
    # Rcpp::List sirt_rcpp_D_eigenvalues( Rcpp::NumericMatrix Xr, int D, int maxit, double conv )
    res <- sirt_rcpp_D_eigenvalues( Xr=X, D=D, maxit=maxit, conv=conv)
    return(res)
}
