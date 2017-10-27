## File Name: isop_tests_cpp.R
## File Version: 0.01


##########################################################
# call to Rcpp function
isop_tests_cpp <- function ( dat , dat.resp , weights , jackunits , JJ ){ 
	isop_tests_C( dat=dat,  dat_resp=dat.resp, weights=weights, 
		jackunits=jackunits,  JJ=JJ )
}	
#############################################################
