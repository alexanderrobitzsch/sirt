## File Name: IRT.likelihood.mirt.R
## File Version: 0.284


IRT_likelihood_mirt_singleGroup <- function( object, type, ... )
{
    object0 <- object
    object <- mirt.wrapper.posterior(mirt.obj=object)
    ll <- object[[ type ]]
    attr(ll,"theta") <- object$theta.k
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- 1
    attr(ll,"pweights") <- object$pweights
    return(ll)
}


#--- likelihood singleGroup
IRT.likelihood.mirt.singleGroup <- function( object, ... )
{
    ll <- IRT_likelihood_mirt_singleGroup( object=object, type="f.yi.qk", ...)
    return(ll)
}
IRT.likelihood.ConfirmatoryClass <- IRT.likelihood.mirt.singleGroup
IRT.likelihood.ExploratoryClass <- IRT.likelihood.mirt.singleGroup
IRT.likelihood.SingleGroupClass <- IRT.likelihood.mirt.singleGroup


#--- likelihood multipleGroup
IRT_likelihood_mirt_multipleGroup <- function( object, type, ... )
{
    mobj <- object
    groups <- object@Data$groupNames
    G <- length(groups)
    prob.theta <- list()
    ll_list <- list()
    ind_group <- list()
    pweights <- list()
    for (gg in 1:G){
        object <- mirt.wrapper.posterior(mirt.obj=mobj, group=groups[gg])
        if (gg==1){
            theta <- object$theta.k
            TP <- object$TP
        }
        prob.theta[[gg]] <- object$pi.k
        ll_list[[gg]] <- object[[type]]
        ind_group[[gg]] <- object$ind_group
        pweights[[gg]] <- object$pweights
    }
    ll <- matrix(NA, nrow=object$N_orig, ncol=object$TP)
    ll_pw <- rep(NA, object$N_orig)
    for (gg in 1:G){
        ll[ ind_group[[gg]], ] <- ll_list[[gg]]
        ll_pw[ ind_group[[gg]] ] <- pweights[[gg]]
    }
    prob.theta <- matrix( unlist(prob.theta), nrow=TP, ncol=G)
    attr(ll,"theta") <- theta
    attr(ll,"prob.theta") <- prob.theta
    attr(ll,"G") <- G
    attr(ll,"pweights") <- ll_pw
    return(ll)
}


IRT.likelihood.mirt.multipleGroup <- function( object, ... )
{
    ll <- IRT_likelihood_mirt_multipleGroup( object=object, type="f.yi.qk", ... )
    return(ll)
}
IRT.likelihood.MultipleGroupClass <- IRT.likelihood.mirt.multipleGroup

