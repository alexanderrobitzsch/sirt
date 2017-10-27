## File Name: md.pattern.sirt.R
## File Version: 0.06

###########################################
# Function for analyzing response patterns
md.pattern.sirt <- function(dat){
    dat <- as.matrix(dat)
	if ( ncol(dat)>1000 ){
	   stop("Function only works for datasets with fewer than 1000 variables!\n")
	}
    res <- md_pattern_rcpp( dat_=dat )
	rp_unique <- unique(res$unique_resp_patt)
	res$unique_resp_patt <- match( res$unique_resp_patt , rp_unique )	
	res$resp_patt <- match( res$resp_patt , rp_unique )
    res$dat.ordered <- res$dat[ order( res$resp_patt ) , ]
    return(res)
}
#*******************************************			
# calling the Rcpp function
md_pattern_rcpp <- function (dat_){ 
	md_pattern_csource( dat_ )
}			
#*******************************************					
