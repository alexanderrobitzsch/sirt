## File Name: summary.isop.test.R
## File Version: 0.07
## File Last Change: 2018-12-30

####################################################
# summary for ISOP test
summary.isop.test <- function( object, ... )
{
    obji <- object$itemstat
    VV <- ncol(obji)
    cat("*** Test for the W1 Axiom in the ISOP Model **** \n\n")
    for (vv in 2:VV){
        obji[,vv] <- round( obji[,vv],3)
    }
    print(obji)
    cat(paste0("\n-- Statistical inference is based on ", object$JJ,
            " jackknife units.\n"))
}
