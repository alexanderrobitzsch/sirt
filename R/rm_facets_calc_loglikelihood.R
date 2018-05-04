## File Name: rm_facets_calc_loglikelihood.R
## File Version: 0.06

rm_facets_calc_loglikelihood <- function( tau.item, a.rater, Qmatrix, b.item, VV, K, I, TP, a.item, 
            b.rater, item.index, rater.index, theta.k, RR, dat2, dat2.resp, pi.k=NULL, dat2.ind.resp,
            mu=NULL, sigma=NULL, b.rater.center , a.rater.center, a.item.center, a_lower, a_upper )  
{
    #--- adjust pi.k probabilities
    if ( ( ! is.null(mu) ) | ( ! is.null(sigma) ) ){
        pi.k <- sirt_dnorm_discrete( x = theta.k, mean=mu, sd = sigma )                
    }
    #--- center parameters
    a.item <- rm_squeeze(x=a.item, lower=a_lower, upper=a_upper )
    a.rater <- rm_squeeze(x=a.rater, lower=a_lower, upper=a_upper )    
    a.rater <- rm_center_vector( vec=a.rater, center_type=a.rater.center, do_log=TRUE )
    a.item <- rm_center_vector( vec=a.item, center_type=a.item.center, do_log=TRUE )
    b.rater <- rm_center_vector( vec=b.rater, center_type=b.rater.center, do_log=FALSE )
    
    #--- calculate probabilities
    probs <- rm_facets_calcprobs( tau.item=tau.item, b.rater=b.rater, Qmatrix=Qmatrix, b.item=b.item, VV=VV, K=K, 
                        I=I, TP=TP, a.item=a.item, a.rater=a.rater, item.index=item.index, 
                        rater.index=rater.index, theta.k=theta.k, RR=RR ) 
                        
    #--- calculate posterior
    res <- rm_posterior( dat2=dat2, dat2.resp=dat2.resp, TP=TP, pi.k=pi.k, K=K, I=I, probs=probs, dat2.ind.resp=dat2.ind.resp ) 

    #--- output
    res <- list(ll=res$ll, a.item=a.item, pi.k=pi.k, a.rater=a.rater, b.rater=b.rater, a.item=a.item)
    return(res)
}
