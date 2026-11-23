## File Name: mgsem_loglike_suffstat_derivative_parameter.R
## File Version: 0.223

mgsem_loglike_suffstat_derivative_parameter <- function(est, dermoments, suffstat,
        type, i1, i2, h, is_B=FALSE, eps=1e-12, num_approx=FALSE)
{

    # num_approx <- TRUE

    #** derivative of moments with respect to parameter
    res <- mgsem_moments_derivative_parameter( est=est, type=type, i1=i1, i2=i2,
                h=h, is_B=is_B, eps=eps, num_approx=num_approx )
    Mu_der <- res$Mu_der
    Sigma_der <- res$Sigma_der
    Sigma_der_logical <- res$Sigma_der_logical

    #** derivatives
    if (res$calc_Mu){
        der1 <- sum( dermoments$dermean * Mu_der )
    } else {
        der1 <- 0
    }
    N <- suffstat$N
    S1 <- dermoments$S1
    y <- dermoments$y
    S2 <- dermoments$S2
    S3 <- dermoments$S3

    if (res$calc_Sigma){

        t4a <- sirt_rcpp_mgsem_loglike_derivative_parameters( S1=S1, S3=S3,
                        y=y, Sigma_der=Sigma_der, Sigma_der_logical=Sigma_der_logical)

        der2 <- sum( -N/2*t4a )
    } else {
        der2 <- 0
    }

    # output
    res <- der1+der2
    return(res)
}
