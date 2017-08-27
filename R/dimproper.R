## File Name: dimproper.R
## File Version: 0.02
## File Last Change: 2017-01-18 11:02:46

###################################################
# improper density which is constant to 1
dimproper <- function(x){
	N <- length(x)
	dx <- rep(1,N)
	return(dx)
}
###################################################
