## File Name: rm_sdt_calc_loglikelihood.R
## File Version: 0.09

rm_sdt_calc_loglikelihood <- function( c.rater, Qmatrix, tau.item, VV, K, I, TP, a.item, 
			d.rater, item.index, rater.index, theta.k, RR, dat2, dat2.resp, pi.k=NULL, dat2.ind.resp,
			mu=NULL, sigma=NULL, center_log_a=FALSE, a_center_type=2)  
{
	#--- adjust pi.k probabilities
	if ( ( ! is.null(mu) ) | ( ! is.null(sigma) ) ){
		pi.k <- sirt_dnorm_discrete( x = theta.k, mean=mu, sd = sigma )				
	}
	#--- center a.item
	if (center_log_a){
		a.item <- rm_center_vector( vec=a.item, center_type=a_center_type, do_log=TRUE )
	}
	#--- calculate probabilities
	res <- rm_hrm_calcprobs( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item, 
						d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR ) 
	probs <- res$prob.total				
	#--- calculate posterior
	res <- rm_posterior( dat2=dat2, dat2.resp=dat2.resp, TP=TP, pi.k=pi.k, K=K, I=I, probs=probs, dat2.ind.resp=dat2.ind.resp ) 
	#--- output
	res <- list(ll=res$ll, a.item=a.item, pi.k=pi.k)
	return(res)
}
