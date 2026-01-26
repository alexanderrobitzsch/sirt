## File Name: linking_sl_optim_function_R.R
## File Version: 0.05


linking_sl_optim_function_R <- function(NI, NS, dist, aM, bM, theta,
    prob_theta, est_pars, wgtM, a, b, mu, sigma, eps, pow=1)
{
    # logit(p)=a*(th-b)
    # th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)
    val <- 0
    for (ss in 1L:NS){

        T1 <- 0
        T0 <- 0

        for (ii in 1L:NI){

            if (est_pars[ii,ss]){
                p_obs <- stats::plogis( aM[ii,ss] * (theta - bM[ii,ss] ) )
                # a_exp <- a[ii] * sigma[ss]
                # b_exp <- ( b[ii] - mu[ss] ) / sigma[ss]
                p_exp <- stats::plogis( a[ii]*sigma[ss]*theta - a[ii]*(b[ii]-mu[ss]))
                T0 <- T0+wgtM[ii,ss]*p_obs
                T1 <- T1+wgtM[ii,ss]*p_exp

            }  # end pars


        }  # end ii

        val <- val + sum( (T1-T0)^2*prob_theta )

    }  # end item ss

    #--- output
    return(val)
}
