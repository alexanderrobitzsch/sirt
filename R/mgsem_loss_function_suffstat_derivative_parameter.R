## File Name: mgsem_loss_function_suffstat_derivative_parameter.R
## File Version: 0.108

mgsem_loss_function_suffstat_derivative_parameter <- function(est, dermoments, suffstat,
        type, i1, i2, h, is_B=FALSE, num_approx=FALSE, eps=1e-12, vech_rcpp=TRUE,
        output_all=FALSE)
{
    #** derivative of moments with respect to parameter
    res <- mgsem_moments_derivative_parameter( est=est, type=type, i1=i1, i2=i2,
                h=h, is_B=is_B, num_approx=num_approx, eps=eps )
    Mu_der <- res$Mu_der
    Sigma_der <- res$Sigma_der
    Sigma_der_logical <- res$Sigma_der_logical
    t1 <- 0
    t2 <- 0

    #*** means
    if (res$calc_Mu){
        t1 <- sum(Mu_der*dermoments$dermean)
    }
    #*** covariance
    if (res$calc_Sigma){
        if (!vech_rcpp){
            vech_sigma_der <- mgsem_vech(x=Sigma_der)
            vech_sigma_der_logical <- mgsem_vech(x=Sigma_der_logical)
        } else {
            vech_sigma_der <- sirt_rcpp_mgsem_vech_numeric(A=Sigma_der)
            vech_sigma_der_logical <- sirt_rcpp_mgsem_vech_logical(A=Sigma_der_logical)
        }
        # t2 <- sum(vech_sigma_der*dermoments$dercov)
        t2 <- sirt_rcpp_mgsem_sumproduct_logical(y=vech_sigma_der,
                    x=dermoments$dercov, y_logical=vech_sigma_der_logical)
    }

    #*** output
    res <- t1+t2
    return(res)
}
