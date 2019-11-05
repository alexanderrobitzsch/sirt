## File Name: IRT.expectedCounts.mirt.R
## File Version: 0.07

# IRT.expectedCounts.mirt

#- counts single group
IRT.expectedCounts.SingleGroupClass <- function( object, ... )
{
    type <- "exp_counts"
    object <- mirt.wrapper.posterior(mirt.obj=object)
    ll <- object[[type]]
    attr(ll,"theta") <- object$theta.k
    attr(ll,"prob.theta") <- object$pi.k
    attr(ll,"G") <- 1
    attr(ll,"pweights") <- object$pweights
    return(ll)
}


#--- counts multipleGroup
IRT.expectedCounts.MultipleGroupClass <- function( object, ... )
{
    mobj <- object
    groups <- object@Data$groupNames
    G <- length(groups)
    prob.theta <- list()
    ll_list <- list()
    ind_group <- list()
    pweights <- list()
    type <- "exp_counts"
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
    dims <- dim(ll_list[[1]])[1:3]
    ll <- array(NA, dim=c(dims, G))
    ll_pw <- rep(NA, object$N_orig)
    for (gg in 1:G){
        ll[,,,gg] <- ll_list[[gg]]
        ll_pw[ ind_group[[gg]] ] <- pweights[[gg]]
    }
    prob.theta <- matrix( unlist(prob.theta), nrow=TP, ncol=G)
    attr(ll,"theta") <- theta
    attr(ll,"prob.theta") <- prob.theta
    attr(ll,"G") <- G
    attr(ll,"pweights") <- ll_pw
    return(ll)
}
