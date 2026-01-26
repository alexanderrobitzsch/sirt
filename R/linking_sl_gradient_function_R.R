## File Name: linking_sl_gradient_function_R.R
## File Version: 0.06


linking_sl_gradient_function_R <- function(NI, NS, dist, aM, bM, theta,
    prob_theta, est_pars, wgtM, a, b, mu, sigma, eps, pow=1,
    index_a, index_b, index_mu, index_sigma, parnames, NP)
{
    # logit(p)=a*(th-b)
    # th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)

    grad <- rep(NA, NP)
    names(grad) <- parnames

    mu_grad <- rep(0, NS)
    sigma_grad <- rep(0, NS)
    a_grad <- rep(0, NI)
    b_grad <- rep(0, NI)

    val <- 0
    for (ss in 1L:NS){

        T1 <- 0
        T0 <- 0

        mu_temp <- 0
        sigma_temp <- 0
        a_temp <- as.list(0*(1L:NI))
        b_temp <- as.list(0*(1L:NI))

        for (ii in 1L:NI){

            if (est_pars[ii,ss]){
                p_obs <- stats::plogis( aM[ii,ss] * (theta - bM[ii,ss] ) )
                # a_exp <- a[ii] * sigma[ss]
                # b_exp <- ( b[ii] - mu[ss] ) / sigma[ss]
                p_exp <- stats::plogis( a[ii]*sigma[ss]*theta - a[ii]*(b[ii]-mu[ss] ))
                T0 <- T0+wgtM[ii,ss]*p_obs
                T1 <- T1+wgtM[ii,ss]*p_exp

                fac_ii <- wgtM[ii,ss]*p_exp*(1-p_exp)

                mu_temp <- mu_temp - fac_ii*a[ii]
                sigma_temp <- sigma_temp - fac_ii*a[ii]*theta
                a_temp[[ii]] <- -fac_ii*( sigma[ss]*theta - ( b[ii] - mu[ss] ) )
                b_temp[[ii]] <- fac_ii*a[ii]

            }  # end pars

        }  # end ii

        mu_grad[ss] <- -sum( (T1-T0)*2*mu_temp*prob_theta)
        sigma_grad[ss] <- -sum( (T1-T0)*2*sigma_temp*prob_theta)

        for (ii in 1L:NI){
            a_grad[ii] <- a_grad[ii] -sum( (T1-T0)*2*a_temp[[ii]]*prob_theta)
            b_grad[ii] <- b_grad[ii] -sum( (T1-T0)*2*b_temp[[ii]]*prob_theta)

        }

        val <- val + sum( (T1-T0)^2*prob_theta )

    }  # end study ss


    grad[ index_mu ] <- mu_grad[-1]
    grad[ index_sigma ] <- sigma_grad[-1]
    grad[ index_a ] <- a_grad
    grad[ index_b ] <- b_grad


    #--- output
    return(grad)
}
