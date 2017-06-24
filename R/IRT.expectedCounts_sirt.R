

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
