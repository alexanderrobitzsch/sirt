## File Name: gom_em_loglike_grad.R
## File Version: 0.278
## File Last Change: 2023-03-15



gom_em_loglike_grad <- function(x, ind_lambda, ind_pi, I, K, theta.k, theta0.k,
        dat2, dat2.resp, TP, model="GOM", h=1e-4, ind_mu, ind_sigma, theta_grid,
        weights, lambda_partable, lambda_index_blocks)
{
    ncat <- 2
    #-- function value
    res <- gom_em_loglike_parameter_conversion( x=x, ind_lambda=ind_lambda,
                ind_pi=ind_pi, I=I, K=K, ind_mu=ind_mu, ind_sigma=ind_sigma, model=model,
                theta_grid=theta_grid, lambda_partable=lambda_partable )
    lambda <- res$lambda
    pi.k <- res$pi.k
    res <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=NULL, theta0.k=theta0.k )
    probs00 <- probs <- res$probs
    probs0 <- sirt_probs_dichotomous_to_array(probs=res$probs)
    fyiqk <- sirt_rcpp_gom_em_likelihood( probs=probs0, ncat=ncat, TP=TP,
                        dat2=dat2, dat2resp=dat2.resp)
    ll <- sirt_rcpp_gom_em_log_likelihood( fyiqk=fyiqk, pik=pi.k, weights=weights)

    #-- gradient with respect to item parameters
    NP <- length(x)
    grad1 <- rep(NA, NP)
    item_vec <- 1:I
    B <- length(lambda_index_blocks)
    for (jj in 1:B){
        ind <- lambda_index_blocks[[jj]]
        lbb <- lambda_partable[ lambda_partable$block==jj,, drop=FALSE]
        items <- lbb[ lbb$free==1, 'item'] - 1
        x1 <- sirt_add_increment(x=x, pos=ind, value=h)
        lambda <- gom_em_loglike_parameter_conversion( x=x1, ind_lambda=ind_lambda,
                        ind_pi=ind_pi, I=I, K=K, ind_mu=ind_mu, ind_sigma=ind_sigma,
                        model=model, theta_grid=theta_grid,
                        lambda_partable=lambda_partable )$lambda
        probs <- gom_em_calc_probs( lambda=lambda, theta.k=theta.k, b=NULL,
                        theta0.k=theta0.k )$probs
        probs_h <- sirt_probs_dichotomous_to_array(probs=probs)
        # update log-likelihood in gradient
        ll2 <- sirt_rcpp_gom_em_loglike_gradient( probs=probs0, probs_h=probs_h,
                    ncat=ncat, TP=TP, dat2=dat2, dat2resp=dat2.resp, pik=pi.k,
                    items=items, fyiqk=fyiqk, weights=weights)
        grad1[ind] <- (ll2 - ll) / h
    }

    #-- gradient with respect to trait distribution parameters
    ind_pi1 <- ind_pi
    if (model=='GOMnormal'){
        ind_pi1 <- c(ind_mu, ind_sigma)
    }
    for (pp in ind_pi1){
        x1 <- sirt_add_increment(x=x, pos=pp, value=h)
        res <- gom_em_loglike_parameter_conversion( x=x1, ind_lambda=ind_lambda,
                    ind_pi=ind_pi, I=I, K=K, ind_mu=ind_mu, ind_sigma=ind_sigma,
                    model=model, theta_grid=theta_grid, lambda_partable=lambda_partable)
        pi.k <- res$pi.k
        ll2 <- sirt_rcpp_gom_em_log_likelihood( fyiqk=fyiqk, pik=pi.k, weights=weights)
        grad1[pp] <- (ll2 - ll) / h
    }
    return(grad1)
}
