## File Name: IRT.expectedCounts_sirt.R
## File Version: 0.02
## File Last Change: 2017-01-18 11:02:47


###########################################################
# object of class xxirt
IRT.expectedCounts.xxirt <- function( object , ... ){    
	ll <- object$n.ik
    attr(ll,"theta") <- object$Theta
	attr(ll,"prob.theta") <- object$probs_Theta
	attr(ll,"G") <- object$G
    return(ll)
        }
###########################################################
