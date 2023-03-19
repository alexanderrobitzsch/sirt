## File Name: linking_haebara_optim_function_R.R
## File Version: 0.171


linking_haebara_optim_function_R <- function(NI, NS, dist, aM, bM, theta,
    prob_theta, est_pars, wgtM, a, b, mu, sigma, eps, pow=1)
{
    # logit(p)=a*(th-b)
    # th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)
    val <- 0
    for (ii in 1:NI){
        for (ss in 1:NS){
            if (est_pars[ii,ss]){
                p_obs <- stats::plogis( aM[ii,ss] * (theta - bM[ii,ss] ) )
                a_exp <- a[ii] * sigma[ss]
                b_exp <- ( b[ii] - mu[ss] ) / sigma[ss]
                # p_exp <- stats::plogis( a_exp * (theta - b_exp ) )
                p_exp <- stats::plogis( a_exp*theta - a[ii]*( b[ii] - mu[ss] ) )
                ## a_exp[ii]=a[ii]*sigma[ss]
                ## a_exp[ii] * b_exp[ii]=a[ii]*(b[ii]-mu[ss])
                dist2 <- (p_obs - p_exp)^2
                if (dist=='L2'){
                    dist1 <- sum( dist2*prob_theta )
                }
                if (dist=='L1'){
                    dist1 <- sum( sqrt( dist2 + eps )*prob_theta )
                }
                if (dist=='Lp'){
                    dist1 <- sum( ( dist2 + eps )^(pow/2)*prob_theta )
                }
                val <- val + wgtM[ii,ss]*dist1
            }
        }
    }
    #--- output
    return(val)
}
