## File Name: rm_grouped_expected_likelihood.R
## File Version: 0.01
## File Last Change: 2017-10-02 12:01:00

rm_grouped_expected_likelihood <- function(pjk, n.ik, diffindex=NULL, eps=1E-80)
{
    ll0 <- rowSums( n.ik * log(pjk+eps) )
	if ( ! is.null(diffindex) ){
		ll0 <- rowsum(ll0, diffindex )[,1]
	}
	return(ll0)
}
