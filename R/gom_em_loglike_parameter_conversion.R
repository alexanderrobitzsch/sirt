## File Name: gom_em_loglike_parameter_conversion.R
## File Version: 0.192


gom_em_loglike_parameter_conversion <- function(x, ind_lambda, ind_pi, I, K,
            ind_mu, ind_sigma, model, theta_grid, lambda_partable)
{
    mu <- NULL
    Sigma <- NULL
    #- item response probabilities
    lambda_logit <- x[ind_lambda]
    lambda_logit <- lambda_logit[ lambda_partable$par_index ]
    lambda <- gom_em_extract_lambda_matrix(lambda_logit=lambda_logit, I=I, K=K)
    # trait distribution
    if ( model=="GOM" ){
        pi_k_logit <- x[ind_pi]
        pi.k <- sirt_logit_to_probs(y=pi_k_logit)
    }
    if (model=="GOMnormal"){
        mu <- x[ind_mu]
        sigma_chol <- x[ind_sigma]
        K1 <- K-1
        sc <- matrix(0, nrow=K1, ncol=K1)
        sc[!upper.tri(sc)] <- sigma_chol
        Sigma <- tcrossprod(x=sc)
        pi.k <- sirt_dmvnorm_discrete(x=theta_grid, mean=mu, sigma=Sigma, eps=1e-10)
    }
    #--- output
    res <- list(lambda=lambda, pi.k=pi.k, mu=mu, Sigma=Sigma)
    return(res)
}
