## File Name: gom_em_loglike_parameter_conversion.R
## File Version: 0.02


gom_em_loglike_parameter_conversion <- function(x, ind_lambda, ind_pi, I, K)
{
    pi_k_logit <- x[ind_pi]
    lambda_logit <- x[ind_lambda]
    pi.k <- sirt_logit_to_probs(y=pi_k_logit)
    lambda <- gom_em_extract_lambda_matrix(lambda_logit=lambda_logit, I=I, K=K)
    #--- output
    res <- list(lambda=lambda, pi.k=pi.k)
    return(res)
}
