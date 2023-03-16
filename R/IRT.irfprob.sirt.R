## File Name: IRT.irfprob.sirt.R
## File Version: 0.16
## File Last Change: 2018-12-30

########################################################
# convert theta.k into a matrix
vec2mat.sirt <- function( theta.k){
    if ( ! is.matrix( theta.k) ){
        theta.k <- matrix(theta.k, ncol=1 )
    }
    return(theta.k)
}
########################################################

########################################################
# irfprob rasch.mml
IRT.irfprob.rasch.mml <- function( object, ... )
{
    ll <- object$rprobs
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- object$G
    return(ll)
}
########################################################

########################################################
# irfprob smirt
IRT.irfprob.smirt <- function( object, ... )
{
    ll <- aperm( object$probs, c(3,2,1) )
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- object$G
    return(ll)
}
########################################################

########################################################
# irfprob rasch.mirtlc
IRT.irfprob.rasch.mirtlc <- function( object, ... )
{
    ll <- object$rprobs
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- object$G
    return(ll)
}
########################################################


########################################################
# irfprob gom
IRT.irfprob.gom <- function( object, ... )
{
    ll <- object$probs
    attr(ll,"theta") <- vec2mat.sirt( object$theta.k )
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- object$G
    return(ll)
}
########################################################
# irfprob rm.facets
IRT.irfprob.rm.facets <- IRT.irfprob.gom
########################################################
# irfprob rm.sdt
IRT.irfprob.rm.sdt <- IRT.irfprob.gom
########################################################
# irfprob prob.guttman
IRT.irfprob.prob.guttman <- IRT.irfprob.gom
########################################################
# irfprob xxirt
IRT.irfprob.xxirt <- function( object, ... )
{
    ll <- object$probs_items
    attr(ll,"theta") <- object$Theta
    attr(ll,"prob.theta") <- object$probs_Theta
    attr(ll,"G") <- object$G
    return(ll)
}
########################################################

