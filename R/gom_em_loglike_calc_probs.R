## File Name: gom_em_loglike_calc_probs.R
## File Version: 0.03
## File Last Change: 2019-05-15

gom_em_loglike_calc_probs <- function(x, ind_pi, ind_lambda, I, K, theta.k,
    theta0.k)
{
    pi_k_logit <- x[ind_pi]
    lambda_logit <- x[ind_lambda]
    pi.k <- sirt_logit_to_probs(y=pi_k_logit)
    lambda <- gom_em_extract_lambda_matrix(lambda_logit=lambda_logit, I=I, K=K)
    res <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=NULL,
                                theta0.k=theta0.k )
    probs <- res$probs
    return(probs)
}
