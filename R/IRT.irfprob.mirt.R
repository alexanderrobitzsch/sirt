## File Name: IRT.irfprob.mirt.R
## File Version: 0.27


#--- irfprob singleGroup
IRT.irfprob.mirt.singleGroup <- function( object, ... )
{
    object <- mirt.wrapper.posterior(mirt.obj=object)
    ll <- object$probs
    attr(ll,"theta") <- object$theta.k
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- 1
    return(ll)
}
IRT.irfprob.ConfirmatoryClass <- IRT.irfprob.mirt.singleGroup
IRT.irfprob.ExploratoryClass <- IRT.irfprob.mirt.singleGroup
IRT.irfprob.SingleGroupClass <- IRT.irfprob.mirt.singleGroup


#--- irfprob multipleGroup
#--- by default, only item response functions from the first group are extracted
IRT.irfprob.mirt.multipleGroup <- function( object, group=1, ... )
{
    mobj <- object
    groups <- object@Data$groupNames
    G <- length(groups)
    prob.theta <- list()
    for (gg in 1:G){
        object <- mirt.wrapper.posterior(mirt.obj=mobj, group=groups[gg])
        if (gg==group){
            ll <- object$probs
            theta <- object$theta.k
            TP <- object$TP
        }
        prob.theta[[gg]] <- object$pi.k
    }
    prob.theta <- matrix( unlist(prob.theta), nrow=TP, ncol=G)
    attr(ll,"theta") <- theta
    attr(ll,"prob.theta") <- prob.theta
    attr(ll,"G") <- G
    return(ll)
}
IRT.irfprob.MultipleGroupClass <- IRT.irfprob.mirt.multipleGroup
