## File Name: rm_calclike.R
## File Version: 0.03
## File Last Change: 2017-10-02 18:53:47

################################################
# calculation of the likelihood
rm_calclike <- function (dat2, dat2resp, probs,K){ 
	RM_CALCPOST( DAT2=dat2, DAT2RESP=dat2resp, PROBS=probs, KK=K)
}
