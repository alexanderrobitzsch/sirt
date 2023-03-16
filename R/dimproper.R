## File Name: dimproper.R
## File Version: 0.04
## File Last Change: 2018-12-30

###################################################
# improper density which is constant to 1
dimproper <- function(x){
    N <- length(x)
    dx <- rep(1,N)
    return(dx)
}
###################################################
