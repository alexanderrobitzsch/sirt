## File Name: gom_em_loglike_opt_fun.R
## File Version: 0.271
## File Last Change: 2022-05-16



gom_em_loglike_opt_fun <- function(x, ind_lambda, ind_pi, I, K, theta.k, theta0.k,
        dat2, dat2.resp, TP, model, ind_mu, ind_sigma, theta_grid, weights,
        lambda_partable, lambda_index_blocks)
{
    ncat <- 2
    res <- gom_em_loglike_parameter_conversion( x=x, ind_lambda=ind_lambda,
                ind_pi=ind_pi, I=I, K=K, ind_mu=ind_mu, ind_sigma=ind_sigma, model=model,
                theta_grid=theta_grid, lambda_partable=lambda_partable )
    lambda <- res$lambda
    pi.k <- res$pi.k
    mu <- res$mu
    Sigma <- res$Sigma

    res <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=NULL, theta0.k=theta0.k )
    #- convert probabilities to array
    probs1 <- sirt_probs_dichotomous_to_array(probs=res$probs)
    #- compute individual likelihood
    fyiqk <- sirt_rcpp_gom_em_likelihood( probs=probs1, ncat=ncat, TP=TP,
                        dat2=dat2, dat2resp=dat2.resp)
    #- compute total log-likelihood
    ll <- sirt_rcpp_gom_em_log_likelihood( fyiqk=fyiqk, pik=pi.k, weights=weights)
    return(ll)
}

