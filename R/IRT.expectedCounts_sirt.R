## File Name: IRT.expectedCounts_sirt.R
## File Version: 0.171



#*** object of class xxirt
IRT.expectedCounts.xxirt <- function( object, ... )
{
    ll <- object$n.ik
    attr(ll,'theta') <- object$Theta
    attr(ll,'prob.theta') <- object$probs_Theta
    attr(ll,'G') <- object$G
    return(ll)
}

#*** object of class rasch.mml
IRT.expectedCounts.rasch.mml <- function( object, ... )
{
    njk <- object$n.jk
    rjk <- object$r.jk
    dims <- dim(njk)
    ll <- array(0, dim=c(dims[1],2,dims[2],dims[3]))
    ll[,2,,] <- rjk
    ll[,1,,] <- njk-rjk
    attr(ll,'theta') <- object$theta.k
    attr(ll,'prob.theta') <- object$pi.k
    attr(ll,'G') <- object$G
    attr(ll,'dimnames') <- list()
    attr(ll,'dimnames')[[1]] <- colnames(object$dat)
    return(ll)
}
