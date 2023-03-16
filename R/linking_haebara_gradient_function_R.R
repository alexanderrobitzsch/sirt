## File Name: linking_haebara_gradient_function_R.R
## File Version: 0.292


linking_haebara_gradient_function_R <- function(NI, NS, dist, aM, bM, theta,
    prob_theta, est_pars, wgtM, a, b, mu, sigma, eps, index_a, index_b,
    index_mu, index_sigma, parnames, NP, pow=1 )
{

    # logit(p)=a*(th-b)
    # th=SIG*TH+MU=> logit(p)=a*(SIG*TH+MU-b)=a*SIG*(TH-(-MU)/SIG-b/SIG)
    grad <- rep(0, 2*NI+2*(NS-1) )
    names(grad) <- parnames
    for (ii in 1:NI){
        for (ss in 1:NS){
            if (est_pars[ii,ss]){
                p_obs <- stats::plogis( aM[ii,ss] * (theta - bM[ii,ss] ) )
                a_exp <- a[ii] * sigma[ss]
                b_exp <- ( b[ii] - mu[ss] ) / sigma[ss]
                p_exp <- stats::plogis( a_exp * (theta - b_exp ) )
                der <- p_exp*(1-p_exp)
                if (dist=="L2"){
                    der_basis <- -2*(p_obs - p_exp)*prob_theta*der
                }
                if (dist=="L1"){
                    diff2 <- p_obs - p_exp
                    dist2 <- diff2^2
                    der_basis <- -(dist2+eps)^(-.5)*diff2*prob_theta*der
                }
                if (dist=="Lp"){
                    diff2 <- p_obs - p_exp
                    dist2 <- diff2^2
                    der_basis <- -pow*(dist2+eps)^(pow/2-1)*diff2*prob_theta*der
                }
                w_t <- wgtM[ii,ss]
                #- a
                der_t <- sigma[ss]*theta + mu[ss] - b[ii]
                ind <- index_a[ii]
                grad[ind] <- grad[ind] + w_t*sum(der_basis*der_t )
                #- b
                der_t <- -a[ii]
                ind <- index_b[ii]
                grad[ind] <- grad[ind] + w_t*sum(der_basis*der_t )
                #- mu
                if (ss>1){
                    der_t <- a[ii]
                    ind <- index_mu[ss-1]
                    grad[ind] <- grad[ind] + w_t*sum(der_basis*der_t )
                }
                #- sigma
                if (ss>1){
                    der_t <- a[ii]*theta
                    ind <- index_sigma[ss-1]
                    grad[ind] <- grad[ind] + w_t*sum(der_basis*der_t )
                }
            }  # end est_parm[ii,ss]
        }  # end ss
    }  # end ii

    #--- output
    return(grad)
}
